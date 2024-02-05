from Crypto.Cipher import AES
from hashlib import sha256
from Crypto.Util.Padding import *
import os

block_size = 16


def enc_key(m, pwd):
    key = pwd2key(pwd)
    f1 = open('key.txt', 'wb')
    f2 = open('message.txt', 'wb')
    f1.write(key[:32])
    f2.write(pad(m.encode(), block_size))
    f1.close()
    f2.close()
    os.system('./encrypt.exe')
    f1 = open('key.txt', 'wb')
    f1.write(key[32:])
    f1.close()
    os.system('./decrypt.exe')
    f1 = open('key.txt', 'wb')
    f1.write(key[:32])
    f1.close()
    os.system('./encrypt.exe')
    f2 = open('message.txt', 'rb')
    return f2.read()


def dec_key(c, pwd):
    key = pwd2key(pwd)
    f1 = open('key.txt', 'wb')
    f2 = open('message.txt', 'wb')
    f1.write(key[:32])
    f2.write(c)
    f1.close()
    f2.close()
    os.system('./decrypt.exe')
    f1 = open('key.txt', 'wb')
    f1.write(key[32:])
    f1.close()
    os.system('./encrypt.exe')
    f1 = open('key.txt', 'wb')
    f1.write(key[:32])
    f1.close()
    os.system('./decrypt.exe')
    f2 = open('message.txt', 'rb')
    b = f2.read()
    tmp = "]".encode()
    for i in range(len(b) - 1, -1, -1):
        if b[i] == b[i - 1] == tmp[0]:
            b = b[:i + 1]
            break
    return b


def pwd2key(pwd):
    m = sha256(pwd.encode()).hexdigest()
    return m.encode()
