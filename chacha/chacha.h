
void printState(unsigned long long state[]);


void num_to_bytes(unsigned long long num, unsigned char bytes[]) {
    for(int i = 0; i < 4; i++) {
        bytes[i] = (num >> (i*8)) & 0xFF;
    }
}

unsigned long long bytes_to_num(unsigned char bytes[]) {
     unsigned long long num = 0;
    for(int i = 0; i < 4; i++) {
        num = num | (( unsigned long long)(bytes[i]) << (i*8));
    }
    return num;
}

unsigned long long key_to_block(unsigned char keys[], int offset) {
     unsigned long long num = 0;
    for(int i = 0; i < 4; i++) {
        num = num | (( unsigned long long)(keys[i + offset]) << (i*8));
    }
    return num;
}

unsigned long long Add(unsigned long long a, unsigned long long b) {
    return (a+b);
    //return (a+b) & 0xFFFFFFFF;
}

unsigned long long Xor(unsigned long long a, unsigned long long b) {
    return (a^b) & 0xFFFFFFFF;
}

unsigned long long Lrot(unsigned long long a, unsigned int rot_bits) {
    return ((a << rot_bits) | (a >> (32-rot_bits))) & 0xFFFFFFFF;
}

void QuaterRound(unsigned long long *a, unsigned long long *b,
                 unsigned long long *c, unsigned long long *d) {

    *a = Add(*a, *b); *d = Xor(*d, *a); *d = Lrot(*d, 16);
    *c = Add(*c, *d); *b = Xor(*b, *c); *b = Lrot(*b, 12);
    *a = Add(*a, *b); *d = Xor(*d, *a); *d = Lrot(*d, 8);
    *c = Add(*c, *d); *b = Xor(*b, *c); *b = Lrot(*b, 7);
}

void inner_block(unsigned long long state[]) {
    // Column
    QuaterRound(&state[0], &state[4], &state[8], &state[12]);
    QuaterRound(&state[1], &state[5], &state[9], &state[13]);
    QuaterRound(&state[2], &state[6], &state[10], &state[14]);
    QuaterRound(&state[3], &state[7], &state[11], &state[15]);

    // Diagonal
    QuaterRound(&state[0], &state[5], &state[10], &state[15]);
    QuaterRound(&state[1], &state[6], &state[11], &state[12]);
    QuaterRound(&state[2], &state[7], &state[8], &state[13]);
    QuaterRound(&state[3], &state[4], &state[9], &state[14]);

}

void chacha20_block(unsigned long long state[], unsigned char key[],
                    unsigned long long counter, unsigned char nonce[]) {

    state[0] = 0x61707865;
    state[1] = 0x3320646e;
    state[2] = 0x79622d32;
    state[3] = 0x6b206574;

    state[4] = key_to_block(key, 0);
    state[5] = key_to_block(key, 4);
    state[6] = key_to_block(key, 8);
    state[7] = key_to_block(key, 12);
    state[8] = key_to_block(key, 16);
    state[9] = key_to_block(key, 20);
    state[10] = key_to_block(key, 24);
    state[11] = key_to_block(key, 28);

    state[12] = counter;

    state[13] = key_to_block(nonce, 0);
    state[14] = key_to_block(nonce, 4);
    state[15] = key_to_block(nonce, 8);

    unsigned long long working_state[16];
    for(int i=0; i<16; i++) {
        working_state[i] = state[i];
    }

    for(int i=0; i<10; i++) {
        inner_block(working_state);
    }
    for(int i=0; i<16; i++) {
        state[i] = Add(state[i], working_state[i]);
    }
}

void serialize_state(unsigned long long state[], unsigned char serialized_key[]) {
    // state: 4 * 16 -> serialize_state: 1 * 64
    unsigned char tmp[4];
    for(int i = 0; i < 16; i++) {
        num_to_bytes(state[i], tmp);
        for(int j=0; j<4; j++) {
            serialized_key[4*i + j] = tmp[j];
        }
    }
}

void chacha20_encrypt(unsigned char key[], unsigned long long counter, unsigned char nonce[],
                        unsigned char plaintext[], unsigned char ciphertext[],
                        unsigned long long msg_length) {
    unsigned long long state[16];
    unsigned char key_stream[64];

    for(int j = 0; j < msg_length/64; j++) {
        chacha20_block(state, key, counter+j, nonce);
        //printState(state);
        serialize_state(state, key_stream);
        for(int i = 0; i < 64; i++) {
            ciphertext[64*j + i] =
                (unsigned char)(plaintext[64*j + i]) ^ (unsigned char)(key_stream[i]);
        }
    }
    if(msg_length % 64 != 0) {
        int j = msg_length / 64;
        chacha20_block(state, key, counter+j, nonce);
        //printState(state);
        serialize_state(state, key_stream);
        for(int i = 0; i < msg_length % 64; i++) {
            ciphertext[64*j + i] =
                (unsigned char)(plaintext[64*j + i]) ^ (unsigned char)(key_stream[i]);
        }
    }
}

void printState(unsigned long long state[]) {
    for(int i=0; i<16; i++) {
        printf("state[%02d]: %08llx\n", i, state[i]);
    }
    puts("");
}
