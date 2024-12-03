import binascii
from math import ceil

# Initial state vector constants for SM3 hash function
INITIAL_VECTOR = [
    0x7380166f, 0x4914b2b9, 0x172442d7, 0xda8a0600,
    0xa96f30bc, 0x163138aa, 0xe38dee4d, 0xb0fb0e4e,
]

# Tj constants used in compression function
# First 16 values use 0x79cc4519, remaining 48 use 0x7a879d8a
CONSTANTS = [0x79cc4519] * 16 + [0x7a879d8a] * 48

bytes_to_list = lambda data: [i for i in data]

def mixing_function_1(a, b, c, round_idx):
    """
    Boolean function used in the compression function
    For rounds 0-15: FF_j(X,Y,Z) = X ⊕ Y ⊕ Z
    For rounds 16-63: FF_j(X,Y,Z) = (X ∧ Y) ∨ (X ∧ Z) ∨ (Y ∧ Z)
    """
    if round_idx < 16:
        return a ^ b ^ c
    return (a & b) | (a & c) | (b & c)

def mixing_function_2(e, f, g, round_idx):
    """
    Boolean function used in the compression function
    For rounds 0-15: GG_j(X,Y,Z) = X ⊕ Y ⊕ Z
    For rounds 16-63: GG_j(X,Y,Z) = (X ∧ Y) ∨ (¬X ∧ Z)
    """
    if round_idx < 16:
        return e ^ f ^ g
    return (e & f) | (~e & g)

def transform_message_word(input_val):
    """
    P1 permutation function
    Performs cyclic shifts and XOR operations on input value
    P1(X) = X ⊕ (X ≪ 9) ⊕ (X ≪ 17)
    """
    rotated_9 = ((input_val << 9) & 0xffffffff) | ((input_val >> 23) & 0xffffffff)
    rotated_17 = ((input_val << 17) & 0xffffffff) | ((input_val >> 15) & 0xffffffff)
    return input_val ^ rotated_9 ^ rotated_17

def process_message_block(state_vector, message_block):
    """
    Main compression function for SM3
    Takes current state vector and a 512-bit message block
    Returns new state vector after processing
    """
    # Convert bytes to 32-bit words
    words = []
    scaling = 0x1000000
    
    # First 16 words are directly from message block
    for i in range(16):
        accumulated = 0
        for k in range(i*4, (i+1)*4):
            accumulated += message_block[k] * scaling
            scaling //= 0x100
        words.append(accumulated)
        scaling = 0x1000000

    # Message expansion - Generate additional 52 words
    for j in range(16, 68):
        term_1 = words[j-16] ^ words[j-9]
        
        # Perform rotations and transformations
        rotated_3 = ((words[j-3] << 15) & 0xffffffff) | ((words[j-3] >> 17) & 0xffffffff)
        rotated_13 = ((words[j-13] << 7) & 0xffffffff) | ((words[j-13] >> 25) & 0xffffffff)
        
        base_val = term_1 ^ rotated_3
        transformed = base_val ^ ((base_val << 15) & 0xffffffff | (base_val >> 17) & 0xffffffff) ^ \
                     ((base_val << 23) & 0xffffffff | (base_val >> 9) & 0xffffffff)
        
        words.append(transformed ^ rotated_13 ^ words[j-6])

    # Generate W' array
    expanded = [words[j] ^ words[j+4] for j in range(64)]
    
    # Initialize working variables with current state
    a, b, c, d, e, f, g, h = state_vector
    
    # Main compression loop - 64 rounds
    for idx in range(64):
        # Calculate intermediate values
        rotated_a = ((a << 12) & 0xffffffff) | ((a >> 20) & 0xffffffff)
        rotated_constant = ((CONSTANTS[idx] << (idx % 32)) & 0xffffffff) | \
                          ((CONSTANTS[idx] >> (32 - idx % 32)) & 0xffffffff)
        
        intermediate = (rotated_a + e + rotated_constant) & 0xffffffff
        ss1 = ((intermediate << 7) & 0xffffffff) | ((intermediate >> 25) & 0xffffffff)
        ss2 = ss1 ^ rotated_a
        
        # Calculate TT1 and TT2
        tt1 = (mixing_function_1(a, b, c, idx) + d + ss2 + expanded[idx]) & 0xffffffff
        tt2 = (mixing_function_2(e, f, g, idx) + h + ss1 + words[idx]) & 0xffffffff
        
        # Perform rotations
        rotated_b = ((b << 9) & 0xffffffff) | ((b >> 23) & 0xffffffff)
        rotated_f = ((f << 19) & 0xffffffff) | ((f >> 13) & 0xffffffff)
        
        # Update working variables
        d = c
        c = rotated_b
        b = a
        a = tt1
        h = g
        g = rotated_f
        f = e
        
        # Transform TT2 for e
        transformed_tt2 = tt2 ^ ((tt2 << 9) & 0xffffffff | (tt2 >> 23) & 0xffffffff) ^ \
                         ((tt2 << 17) & 0xffffffff | (tt2 >> 15) & 0xffffffff)
        e = transformed_tt2
        # Ensure all values are 32-bit
        a, b, c, d, e, f, g, h = map(lambda x: x & 0xFFFFFFFF, [a, b, c, d, e, f, g, h])
    
    # Combine new state with previous state
    result_state = [a, b, c, d, e, f, g, h]
    return [result_state[i] ^ state_vector[i] for i in range(8)]

def sm3_hash(input_data):
    """
    Main SM3 hash function
    Takes input data as bytes and returns hash value as hex string
    """
    
    data_length = len(input_data)
    remainder = data_length % 64
    message = input_data[:]
    # Append padding bits
    message.append(0x80)
    # Calculate padding length
    padding_length = 56 - remainder - 1 if remainder < 56 else 120 - remainder - 1
    message.extend([0] * padding_length)
    
    # Append message length in bits
    bit_length = data_length * 8
    length_bytes = []
    for _ in range(8):
        length_bytes.insert(0, bit_length & 0xFF)
        bit_length >>= 8
    
    message.extend(length_bytes)
    block_count = len(message) // 64
    
    # Process message in 512-bit blocks
    blocks = [message[i:i+64] for i in range(0, len(message), 64)]
    current_state = INITIAL_VECTOR
    
    for block in blocks:
        current_state = process_message_block(current_state, block)
    
    return ''.join(f'{x:08x}' for x in current_state)

def sm3_kdf(input_str, key_length):
    """
    Key Derivation Function (KDF) based on SM3 hash
    Used to generate keys of specified length from input string
    """
    
    counter = 1
    iterations = ceil(key_length/32)
    input_bytes = [i for i in bytes.fromhex(input_str.decode('utf8'))]
    accumulated = ""
    
    # Generate required number of hash values
    for _ in range(iterations):
        counter_bytes = [i for i in binascii.a2b_hex(f'{counter:08x}'.encode('utf8'))]
        accumulated += sm3_hash(input_bytes + counter_bytes)
        counter += 1
    
    return accumulated[:key_length * 2]