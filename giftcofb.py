                                        # ENCRYPTION GIFT_COFB #
import numpy
from giftbox import *

pad_seq = bitarray([0] * 64)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pad function


def pad(s):
    m = s
    if len(s) != 0 and len(s) % 128 == 0:
        return s

    else:
        m = m + "1"
        for num in range(0, 128 - len(s) % 128 ):
            m = m + "0"

    return m


def rand_key(p):
    key1 = ""

    for i in range(p):
        temp = str(random.randint(0, 1))
        key1 += temp

    return key1


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# G function


def G(str):
    left = str[0:64]
    right = str[64:128]
    left = leftrotate(left, 1)
    state = right + left
    return state


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Encrypt-GIFT-COFB


def encrypt_gift_cofb(n, a, m, k):
    # Pad associated data
    a_temp = pad(a)

    # divide the associated data into blocks of 128 bits
    a_size = int(len(a_temp) / 128)
    m_size= 0

    A = [bitarray] * (a_size)

    temp = 0
    for num in range(0, a_size):
        A[num] = a_temp[temp:temp + 128]
        temp += 128

    # Pad message
    if len(m) != 0:
        m_temp = pad(m)
        m_size = int(len(m_temp)/128)
        # divide the message into blocks of 128 bits
        M = [bitarray] * (m_size)
        temp = 0
        for num in range(0, m_size):
            M[num] = m_temp[temp:temp + 128]
            temp += 128

    X = [bitarray()] * (a_size + m_size + 1)
    Y = [bitarray()] * (a_size + m_size + 1)
    C = [bitarray()] * (m_size + 1)


    Y[0] = GIFT_box(n, k)
    L = Y[0][0:64]


    for i in range(1, a_size):


        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_2, numpy.poly1d(L.tolist())), primitive_pol)


        r = []
        for f in range(0,len(remainder)):
            r.append(int(remainder[f]))


        for j in range(0, len(r)):
             r[j] = r[j] % 2

        if len(r) < 64 :
            for f in range(0, 64 - len(r)):
                r.insert(0,0)

        L = bitarray(r)

        X[i] = A[i - 1] ^ G(Y[i - 1]) ^ (L + pad_seq)
        Y[i] = GIFT_box(X[i], k)




    if len(A) % 128 == 0 and len(A) != 0:
        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())), primitive_pol)

        r = []
        for f in range(0,len(remainder)):
            r.append(int(remainder[f]))


        for j in range(0, len(r)):
             r[j] = r[j] % 2

        if len(r) < 64 :
            for f in range(0, 64 - len(r)):
                r.insert(0,0)

        L = bitarray(r)



    else:
        for n in range(0, 2):
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)

    if(len(m) == 0):
        for n in range(0, 2):
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)

    X[a_size] = A[a_size - 1] ^ G(Y[a_size - 1]) ^ (L + pad_seq)
    Y[a_size] = GIFT_box(X[a_size], k)

    for i in range(1, m_size):

        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_2, numpy.poly1d(L.tolist())),
                                            primitive_pol)

        r = []
        for f in range(0, len(remainder)):
            r.append(int(remainder[f]))

        for j in range(0, len(r)):
            r[j] = r[j] % 2

        if len(r) < 64:
            for f in range(0, 64 - len(r)):
                r.insert(0, 0)

        L = bitarray(r)

        C[i] = M[i - 1] ^ Y[i + a_size - 1]
        X[i + a_size] = M[i - 1] ^ G(Y[i + a_size - 1]) ^ (L + pad_seq)
        Y[i + a_size] = GIFT_box(X[i + a_size], k)


    if len(m) != 0:
        if len(m) % 128 == 0:
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)


        else:

            for n in range (0, 2):
                quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                        primitive_pol)

                r = []
                for f in range(0, len(remainder)):
                    r.append(int(remainder[f]))

                for j in range(0, len(r)):
                    r[j] = r[j] % 2

                if len(r) < 64:
                    for f in range(0, 64 - len(r)):
                        r.insert(0, 0)

                L = bitarray(r)

        C[m_size] = M[m_size - 1] ^ Y[a_size + m_size - 1]
        X[a_size + m_size] = M[m_size - 1] ^ G(Y[a_size + m_size - 1] ^ (L + pad_seq))
        Y[m_size + a_size] = GIFT_box(X[m_size + a_size], k)
        c = bitarray()
        t = bitarray()

        # retrieve cipher text
        for h in range(1, m_size + 1):
            c += C[h]


        c = c[0:len(m)]


        # retrieve tag
        t = Y[a_size + m_size][0:128]

    else:
        c = bitarray()
        t = bitarray()
        t = Y[a_size][0:128]

    list = []
    list.append(c)
    list.append(t)
    return list


                                            # DECRYPTION GIFT_COFB #

def decrypt_gift_cofb(n, a, c, t, k):
    # Pad associated data
    a_temp = pad(a)

    # divide the associated data into blocks of 128 bits
    a_size = int(len(a_temp) / 128)
    c_size = 0

    A = [bitarray] * (a_size)

    temp = 0
    for num in range(0, a_size):
        A[num] = a_temp[temp:temp + 128]
        temp += 128

    # Pad message
    if len(c) != 0:
        c_temp = pad(c)
        c_size = int(len(c_temp) / 128)
        # divide the message into blocks of 128 bits
        C = [bitarray] * (c_size)
        temp = 0
        for num in range(0, c_size):
            C[num] = c_temp[temp:temp + 128]
            temp += 128

    X = [bitarray()] * (a_size + c_size + 1)
    Y = [bitarray()] * (a_size + c_size + 1)
    M = [bitarray()] * (c_size + 1)

    Y[0] = GIFT_box(n, k)
    L = Y[0][0:64]

    for i in range(1, a_size):

        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_2, numpy.poly1d(L.tolist())),
                                            primitive_pol)

        r = []
        for f in range(0, len(remainder)):
            r.append(int(remainder[f]))

        for j in range(0, len(r)):
            r[j] = r[j] % 2

        if len(r) < 64:
            for f in range(0, 64 - len(r)):
                r.insert(0, 0)

        L = bitarray(r)
        X[i] = A[i - 1] ^ G(Y[i - 1]) ^ (L + pad_seq)
        Y[i] = GIFT_box(X[i], k)

    if len(A) % 128 == 0 and len(A) != 0:
        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                            primitive_pol)

        r = []
        for f in range(0, len(remainder)):
            r.append(int(remainder[f]))

        for j in range(0, len(r)):
            r[j] = r[j] % 2

        if len(r) < 64:
            for f in range(0, 64 - len(r)):
                r.insert(0, 0)

        L = bitarray(r)



    else:
        for n in range(0, 2):
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)

    if (len(c) == 0):
        for n in range(0, 2):
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)

    X[a_size] = A[a_size - 1] ^ G(Y[a_size - 1]) ^ (L + pad_seq)
    Y[a_size] = GIFT_box(X[a_size], k)


    for i in range(1, c_size):

        quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_2, numpy.poly1d(L.tolist())),
                                            primitive_pol)

        r = []
        for f in range(0, len(remainder)):
            r.append(int(remainder[f]))

        for j in range(0, len(r)):
            r[j] = r[j] % 2

        if len(r) < 64:
            for f in range(0, 64 - len(r)):
                r.insert(0, 0)

        L = bitarray(r)


        M[i] = C[i - 1] ^ Y[i + a_size - 1]
        X[i + a_size] = M[i] ^ G(Y[i + a_size - 1]) ^ (L + pad_seq)
        Y[i + a_size] = GIFT_box(X[i + a_size], k)

    if len(c) != 0:
        if len(c) % 128 == 0:
            quotient, remainder = numpy.polydiv(numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                                                primitive_pol)

            r = []
            for f in range(0, len(remainder)):
                r.append(int(remainder[f]))

            for j in range(0, len(r)):
                r[j] = r[j] % 2

            if len(r) < 64:
                for f in range(0, 64 - len(r)):
                    r.insert(0, 0)

            L = bitarray(r)
            M[c_size] = C[c_size - 1] ^ Y[a_size + c_size - 1]

        else:
            for n in range(0, 2):
                quotient, remainder = numpy.polydiv(
                    numpy.polynomial.polynomial.polymul(pol_3, numpy.poly1d(L.tolist())),
                    primitive_pol)

                r = []
                for f in range(0, len(remainder)):
                    r.append(int(remainder[f]))

                for j in range(0, len(r)):
                    r[j] = r[j] % 2

                if len(r) < 64:
                    for f in range(0, 64 - len(r)):
                        r.insert(0, 0)

                L = bitarray(r)

            c_dash = len(c) % 128

            sequence = bitarray()
            sequence.append(1)
            for num in range(1, 128 - c_dash):
                sequence.append(0)
            M[c_size] = (C[c_size - 1] ^ Y[a_size + c_size - 1])[0: c_dash] + sequence

        X[a_size + c_size] = M[c_size] ^ G(Y[a_size + c_size - 1] ^ (L + pad_seq))
        Y[c_size + a_size] = GIFT_box(X[c_size + a_size], k)

        m = bitarray()


        # retrieve plain text
        for h in range(1, c_size + 1):
            m += M[h]

        m = m[0:len(c)]

        t_dash = bitarray()

        # calculate tag
        t_dash = Y[a_size + c_size][0:128]


    else:

        m = bitarray()
        t_dash = bitarray()
        t_dash = Y[a_size][0:128]


    if t_dash == t:
        return m

    else:
        return "Error"

# nonce
n = noncegenerator()
# associated data
a = "1010111011110101011010"
# message
m = "0000110110111010110010101011110111010001111011100111101011111000000101000000010101010111101010011001101110011101110000110111000100111000011000000001000010001100111101001110111011101010000101100111010011011001001111110100100000100000100001010110101000111111"
# key
k = "10010100100101001001101011111000111100110100110010010111001011100011000010110000001111101100101010000100010101110010000100010101"

#Encrypting- returns [ciphertext, tag]
t = encrypt_gift_cofb(bitarray(n), bitarray(a), bitarray(m), bitarray(k))
#Decrypting
g  = decrypt_gift_cofb(bitarray(n), bitarray(a), t[0], t[1], bitarray(k))
#checking if the message recovered is equal to the original message
print(g == bitarray(m))
#passing wrong associated data gives error proving that authenticity is ensured
aa = "1010111011110101011000"
c  = decrypt_gift_cofb(bitarray(n), bitarray(aa), t[0], t[1], bitarray(k))
print(c)