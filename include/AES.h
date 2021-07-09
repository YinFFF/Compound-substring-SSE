#ifndef AES_H
#define AES_H

int AES_encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key,
            unsigned char *iv, unsigned char *ciphertext);

int AES_decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
            unsigned char *iv, unsigned char *plaintext);

int HmacEncode(const char * algo, 
		const char * key, unsigned int key_length, 
		const char * input, unsigned int input_length, 
		unsigned char * &output, unsigned int &output_length);


#endif
