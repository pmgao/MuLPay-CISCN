import binascii
from math import ceil

INITIAL_VECTOR = [
    0x7380166f, 0x4914b2b9, 0x172442d7, 0xda8a0600,
    0xa96f30bc, 0x163138aa, 0xe38dee4d, 0xb0fb0e4e,
]

CONSTANTS = [0x79cc4519] * 16 + [0x7a879d8a] * 48

def mixing_function_1(a, b, c, round_idx):
    if round_idx < 16:
        return a ^ b ^ c
    return (a & b) | (a & c) | (b & c)

def mixing_function_2(e, f, g, round_idx):
    if round_idx < 16:
        return e ^ f ^ g
    return (e & f) | (~e & g)

def transform_message_word(input_val):
    rotated_9 = ((input_val << 9) & 0xffffffff) | ((input_val >> 23) & 0xffffffff)
    rotated_17 = ((input_val << 17) & 0xffffffff) | ((input_val >> 15) & 0xffffffff)
    return input_val ^ rotated_9 ^ rotated_17

def process_message_block(state_vector, message_block):
    words = []
    scaling = 0x1000000
    
    for i in range(16):
        accumulated = 0
        for k in range(i*4, (i+1)*4):
            accumulated += message_block[k] * scaling
            scaling //= 0x100
        words.append(accumulated)
        scaling = 0x1000000

    for j in range(16, 68):
        term_1 = words[j-16] ^ words[j-9]
        
        rotated_3 = ((words[j-3] << 15) & 0xffffffff) | ((words[j-3] >> 17) & 0xffffffff)
        rotated_13 = ((words[j-13] << 7) & 0xffffffff) | ((words[j-13] >> 25) & 0xffffffff)
        
        base_val = term_1 ^ rotated_3
        transformed = base_val ^ ((base_val << 15) & 0xffffffff | (base_val >> 17) & 0xffffffff) ^ \
                     ((base_val << 23) & 0xffffffff | (base_val >> 9) & 0xffffffff)
        
        words.append(transformed ^ rotated_13 ^ words[j-6])

    expanded = [words[j] ^ words[j+4] for j in range(64)]
    
    a, b, c, d, e, f, g, h = state_vector
    
    for idx in range(64):
        rotated_a = ((a << 12) & 0xffffffff) | ((a >> 20) & 0xffffffff)
        rotated_constant = ((CONSTANTS[idx] << (idx % 32)) & 0xffffffff) | \
                          ((CONSTANTS[idx] >> (32 - idx % 32)) & 0xffffffff)
        
        intermediate = (rotated_a + e + rotated_constant) & 0xffffffff
        ss1 = ((intermediate << 7) & 0xffffffff) | ((intermediate >> 25) & 0xffffffff)
        ss2 = ss1 ^ rotated_a
        
        tt1 = (mixing_function_1(a, b, c, idx) + d + ss2 + expanded[idx]) & 0xffffffff
        tt2 = (mixing_function_2(e, f, g, idx) + h + ss1 + words[idx]) & 0xffffffff
        
        rotated_b = ((b << 9) & 0xffffffff) | ((b >> 23) & 0xffffffff)
        rotated_f = ((f << 19) & 0xffffffff) | ((f >> 13) & 0xffffffff)
        
        d = c
        c = rotated_b
        b = a
        a = tt1
        h = g
        g = rotated_f
        f = e
        
        transformed_tt2 = tt2 ^ ((tt2 << 9) & 0xffffffff | (tt2 >> 23) & 0xffffffff) ^ \
                         ((tt2 << 17) & 0xffffffff | (tt2 >> 15) & 0xffffffff)
        e = transformed_tt2
        
        a, b, c, d, e, f, g, h = map(lambda x: x & 0xFFFFFFFF, [a, b, c, d, e, f, g, h])
    
    result_state = [a, b, c, d, e, f, g, h]
    return [result_state[i] ^ state_vector[i] for i in range(8)]

def sm3_hash(input_data):
    data_length = len(input_data)
    remainder = data_length % 64
    message = input_data[:]
    message.append(0x80)
    
    padding_length = 56 - remainder - 1 if remainder < 56 else 120 - remainder - 1
    message.extend([0] * padding_length)
    
    bit_length = data_length * 8
    length_bytes = []
    for _ in range(8):
        length_bytes.insert(0, bit_length & 0xFF)
        bit_length >>= 8
    
    message.extend(length_bytes)
    block_count = len(message) // 64
    
    blocks = [message[i:i+64] for i in range(0, len(message), 64)]
    current_state = INITIAL_VECTOR
    
    for block in blocks:
        current_state = process_message_block(current_state, block)
    
    return ''.join(f'{x:08x}' for x in current_state)

def sm3_kdf(input_str, key_length):
    counter = 1
    iterations = ceil(key_length/32)
    input_bytes = [i for i in bytes.fromhex(input_str.decode('utf8'))]
    accumulated = ""
    
    for _ in range(iterations):
        counter_bytes = [i for i in binascii.a2b_hex(f'{counter:08x}'.encode('utf8'))]
        accumulated += sm3_hash(input_bytes + counter_bytes)
        counter += 1
    
    return accumulated[:key_length * 2]
