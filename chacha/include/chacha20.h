#pragma once
#include "array_interface.h"
#include <iostream>

void state_init(uint32_t*, uint8_t*, uint32_t, uint8_t*);
void chacha20_encrypt(uint8_t*, uint8_t*, uint32_t, uint8_t*, uint8_t*, size_t);
//void chacha20_aead_encrypt(uint8_t*, uint8_t*, uint8_t*, size_t, uint8_t*, uint8_t*, uint8_t*, size_t);
//bool chacha20_aead_decrypt(uint8_t*, uint8_t*, uint8_t*, uint8_t*, size_t, uint8_t*, uint8_t*, size_t);
void chacha20_aead_encrypt(uint8_t*, uint8_t*, const array_t*, uint8_t*, uint8_t*, const array_t*);
bool chacha20_aead_decrypt(uint8_t*, uint8_t*, uint8_t*, const array_t*, uint8_t*, const array_t*);

//void chacha20_encrypt(ciphertext, key, counter, nonce, plaintext);

inline void QuaterRound(uint32_t state[16], 
        uint8_t i1, uint8_t i2, uint8_t i3, uint8_t i4) {
    // i1 = a, i2 = b, 
    // i3 = c, i4 = d
    state[i1] += state[i2];
    state[i4] ^= state[i1];
    state[i4] = (state[i4] << 16) | (state[i4] >> 16);

    state[i3] += state[i4];
    state[i2] ^= state[i3];
    state[i2] = (state[i2] << 12) | (state[i2] >> 20);

    state[i1] += state[i2];
    state[i4] ^= state[i1];
    state[i4] = (state[i4] << 8) | (state[i4] >> 24);

    state[i3] += state[i4];
    state[i2] ^= state[i3];
    state[i2] = (state[i2] << 7) | (state[i2] >> 25);
}

