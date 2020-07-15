#pragma once

#include "chacha20.h"
#include "poly1305.h"
#include "array_interface.h"
#include <iostream>
#include <string.h>

// https://tools.ietf.org/html/rfc7539#section-2.5


void state_init(uint32_t state[16], 
        uint8_t key[32], uint32_t count, uint8_t nonce[12]) {
    uint32_t dummy_state[16];
    dummy_state[0] = 0x61707865;
    dummy_state[1] = 0x3320646e;
    dummy_state[2] = 0x79622d32;
    dummy_state[3] = 0x6b206574;
    dummy_state[12] = count;
    memcpy(dummy_state + 4, key, 32);
    memcpy(dummy_state + 13, nonce, 12);

    memcpy(state, dummy_state, 64);
    for (size_t i = 0; i < 10; i++) {
        QuaterRound(dummy_state, 0, 4, 8, 12);
        QuaterRound(dummy_state, 1, 5, 9, 13);
        QuaterRound(dummy_state, 2, 6, 10, 14);
        QuaterRound(dummy_state, 3, 7, 11, 15);
        QuaterRound(dummy_state, 0, 5, 10, 15);
        QuaterRound(dummy_state, 1, 6, 11, 12);
        QuaterRound(dummy_state, 2, 7, 8, 13);
        QuaterRound(dummy_state, 3, 4, 9, 14);
    }

    for (size_t i = 0; i < 16; i++) {
        state[i] += dummy_state[i];
    }
#ifdef NDEBUG
    dump(state, 16);
#endif
}

void chacha20_encrypt(uint8_t *ciphertext, uint8_t key[32], uint32_t counter, 
        uint8_t nonce[12], uint8_t *plaintext, size_t length_plaintext) {
    size_t j;
    uint32_t state[16];

    for (j = 0; j < (length_plaintext + 63) / 64 - 1; j++) {
        state_init(state, key, counter + j, nonce);

        for (size_t k = 0; k < 64; k++) {
            ciphertext[64*j + k] = ((uint8_t*)state)[k] ^ plaintext[64*j + k];
        }
    }

    if (length_plaintext % 64 != 0) {
        //j = (length_plaintext + 63) / 64 - 1;
        state_init(state, key, counter + j, nonce);

        for (size_t k = 0; k < length_plaintext % 64; k++) {
            ciphertext[64*j + k] = ((uint8_t*)state)[k] ^ plaintext[64*j + k];
        }
    }
#ifdef NDEBUG
    puts("");
    dump(ciphertext, length_plaintext);
#endif
}

void poly1305_mac(uint8_t tag[16], uint8_t* msg, size_t length_msg, uint8_t key[32]) {
    mp_limb_t r[3] = {0};
    mp_limb_t s[3] = {0};

    memcpy(s, key + 16, 16);
    memcpy(r, key, 16);
    poly1305aes_test_clamp((uint8_t*)r);

#ifdef NDEBUG
    puts("");
    puts("[*] s:");
    printf("%08lx ", s[1]);
    printf("%08lx\n", s[0]);
    puts("");
    puts("[*] r:");
    printf("%016lx ", r[1]);
    printf("%016lx\n\n", r[0]);
#endif

    mp_limb_t n[YKM_POLY1305_SIZE];
    mp_limb_t acc[3] = {0};

    size_t j;
    for (j = 0; j < length_msg/16; j++) {
        memset(n, 0, sizeof(mp_limb_t)*YKM_POLY1305_SIZE);
        memcpy(n, msg + 16*j, 16);
        n[2] |= 0x01;

        add_3(acc, acc, n);
        DivMod(acc, acc, r);
    }

    if (length_msg % 16 != 0) {
        memset(n, 0, sizeof(mp_limb_t)*YKM_POLY1305_SIZE);
        memcpy(n, msg + 16*j, length_msg % 16);
        ((uint8_t*)n)[length_msg % 16] |= 0x01;

        add_3(acc, acc, n);
        DivMod(acc, acc, r);
    }
    add_3(acc, acc, s);
    memcpy(tag, acc, 16);
}

void poly1305_key_gen(uint8_t poly_key[32], uint8_t key[32], uint8_t nonce[12]) {
    uint32_t state[16];
    state_init(state, key, 0, nonce);
    memcpy(poly_key, state, 32);
}

void chacha20_aead_encrypt(uint8_t* ciphertext, uint8_t* tag, const array_t *aad, 
        uint8_t* key, uint8_t* nonce, const array_t *plaintext) {
    // nonce = constant | iv
    uint8_t otk[32];

    size_t length_mac_data = plaintext->_size + aad->_size;
    length_mac_data += (16 - (plaintext->_size % 16));
    length_mac_data += (16 - (aad->_size % 16));

#if 1
    uint8_t mac_data[length_mac_data + 16];
    memset(mac_data, 0, length_mac_data + 16);
#else
    uint8_t mac_data[length_mac_data + 16] = {0}; 
#endif

    poly1305_key_gen(otk, key, nonce);
    chacha20_encrypt(ciphertext, key, 1, nonce, plaintext->_value, plaintext->_size);

    memcpy(mac_data, aad->_value, aad->_size);
    memcpy(mac_data + (aad->_size + (16 - (aad->_size % 16))), ciphertext, plaintext->_size);
    *(size_t*)(mac_data + length_mac_data) = aad->_size;
    *(size_t*)(mac_data + length_mac_data + 8) = plaintext->_size;

#ifdef NDEBUG
    puts("[*] mac_data");
    for (size_t i = 0; i < length_mac_data + 16; i++) {
        printf("%02x ", mac_data[i]);
        if ((i + 1) % 16 == 0) puts("");
    }
    puts("");
    printf("mac_data length: %ld\n", length_mac_data + 16);
#endif

    poly1305_mac(tag, mac_data, length_mac_data + 16, otk);
}

bool chacha20_aead_decrypt(uint8_t *plaintext, uint8_t *key, 
        uint8_t *nonce, const array_t *aad, uint8_t *tag, const array_t *ciphertext) {
    uint8_t otk[32];
    uint8_t calc_tag[16];

    size_t length_mac_data = ciphertext->_size + aad->_size;
    length_mac_data += (16 - (ciphertext->_size % 16));
    length_mac_data += (16 - (aad->_size % 16));

#if 1
    uint8_t mac_data[length_mac_data + 16];
    memset(mac_data, 0, length_mac_data + 16);
#else
    uint8_t mac_data[length_mac_data + 16] = {0}; 
#endif

    poly1305_key_gen(otk, key, nonce);
    chacha20_encrypt(plaintext, key, 1, nonce, ciphertext->_value, ciphertext->_size);

    memcpy(mac_data, aad->_value, aad->_size);
    memcpy(mac_data + (aad->_size + (16 - (aad->_size % 16))), ciphertext->_value, ciphertext->_size);
    *(size_t*)(mac_data + length_mac_data) = aad->_size; 
    *(size_t*)(mac_data + length_mac_data + 8) = ciphertext->_size;
    poly1305_mac(calc_tag, mac_data, length_mac_data + 16, otk);

    return memcmp(calc_tag, tag, 16) == 0; 
}
