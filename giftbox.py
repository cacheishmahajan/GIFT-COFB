from bitarray import bitarray
import random


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------# ------------------------------------------------------------------------------------------------------------------------------------
# Mechanism for generating nonce


list = []


def noncegenerator():
    while 1:
        nonce = randbit(128)
        if nonce not in list:
            list.append(nonce)
            break

    return nonce


def randbit(n):
    t = bitarray()
    for h in range (0, n):
        t.append(random.getrandbits(1))

    return t

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Polynomials


pol_2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
pol_3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
primitive_pol= [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Permutation Matrix


permutation = [
    [0, 4, 8, 12, 16, 20, 24, 28, 3, 7, 11, 15, 19, 23, 27, 31, 2, 6, 10, 14, 18, 22, 26, 30, 1, 5, 9, 13, 17, 21, 25,
     29],
    [1, 5, 9, 13, 17, 21, 25, 29, 0, 4, 8, 12, 16, 20, 24, 28, 3, 7, 11, 15, 19, 23, 27, 31, 2, 6, 10, 14, 18, 22, 26,
     30],
    [2, 6, 10, 14, 18, 22, 26, 30, 1, 5, 9, 13, 17, 21, 25, 29, 0, 4, 8, 12, 16, 20, 24, 28, 3, 7, 11, 15, 19, 23, 27,
     31],
    [3, 7, 11, 15, 19, 23, 27, 31, 2, 6, 10, 14, 18, 22, 26, 30, 1, 5, 9, 13, 17, 21, 25, 29, 0, 4, 8, 12, 16, 20, 24,
     28]]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Round Constants


round_constants = ["01", "03", "07", "0F", "1F", "3E", "3D", "3B", "37", "2F", "1E", "3C", "39", "33", "27", "0E", "1D",
                   "3A",
                   "35", "2B", "16", "2C", "18", "30", "21", "02", "05", "0B", "17", "2E", "1C", "38", "31", "23", "06",
                   "0D",
                   "1B", "36", "2D", "1A"]


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Encryption box - GIFT cipherllp[[-p


def GIFT_box(plaintext, key):
    cipherstate = initialise_pt(plaintext)
    keystate = initialise_ks(key)

    for i in range(1, 41):
        cipherstate = subcells(cipherstate)
        cipherstate = permBits(cipherstate)
        AddRoundKey(cipherstate, keystate, i)
        keystateupdate(keystate)


    y = finalise(cipherstate)
    return y



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialise Plaintext


def initialise_pt(plaintext):
    # Initialise PlainTextState



     s0 = plaintext[0:32]
     s1 = plaintext[32:64]
     s2 = plaintext[64:96]
     s3 = plaintext[96:128]

     state  = []
     state.append(s0)
     state.append(s1)
     state.append(s2)
     state.append(s3)

     return state


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialise KeyState


def initialise_ks(keystate):
    # Initialise KeyState

    ks0 = keystate[0:16]
    ks1 = keystate[16:32]
    ks2 = keystate[32:48]
    ks3 = keystate[48:64]
    ks4 = keystate[64:80]
    ks5 = keystate[80:96]
    ks6 = keystate[96:112]
    ks7 = keystate[112:128]

    state = []
    state.append(ks0)
    state.append(ks1)
    state.append(ks2)
    state.append(ks3)
    state.append(ks4)
    state.append(ks5)
    state.append(ks6)
    state.append(ks7)

    return state


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SubCells


def subcells(state):


    state[1] = state[1] ^ (state[0] & state[2])
    state[0] = state[0] ^ (state[1] & state[3])
    state[2] = state[2] ^ (state[0] | state[1])
    state[3] = state[3] ^ state[2]
    state[1] = state[1] ^ state[3]
    state[3] = ~state[3]
    state[2] = state[2] ^ (state[0] & state[1])

    list = []
    for i in range (0, 4):
        list.append(state[3 - i])

    return list


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PermBits


def permBits(state):
    tempstate = state
    for count1 in range(0, 4):
        for count2 in range(0, 32):
            tempstate[count1][count2] = state[count1][permutation[count1][count2]]

    return tempstate


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# AddRoundKey


def AddRoundKey(cipherstate, keystate, round):
    # concatenate W6 with W7
    temp1 = keystate[6] + keystate[7]

    # concatenate W2 with W3
    temp2 = keystate[2] + keystate[3]

    cipherstate[1] = cipherstate[1] ^ temp1
    cipherstate[2] = cipherstate[2] ^ temp2

    Const_round = "800000" + round_constants[round - 1]

    # convert hex to binary

    int_value = int(Const_round, base=16)

    binary_value = str(bin(int_value))[2:]


    binary = bitarray(binary_value)
    cipherstate[3] = cipherstate[3] ^ binary

    return


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# KeyStateUpdate


def keystateupdate(state):
    # performing updates
    state[6] = rightrotate(state[6], 2)
    state[7] = rightrotate(state[7], 12)

    state[0] = state[6]
    state[1] = state[7]
    state[2] = state[0]
    state[3] = state[1]
    state[4] = state[2]
    state[5] = state[3]
    state[6] = state[4]
    state[7] = state[5]

    return


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finalise


def finalise(ciphertext):
    list = []

    for count1 in range(0, 4):
        for count2 in range(0, 32):
            list.append(ciphertext[count1][count2])

    return bitarray(list)


# Bit rotation
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



def leftrotate(list, n):
    for h in range(0, n):
        last = list[0]
        for i in range(1, len(list)):
            list[i - 1] = list[i]
        list[len(list) - 1] = last

    return list


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def rightrotate(list, n):
    for h in range(0, n):
        first = list[len(list) - 1]
        for i in range(len(list) - 2, -1, -1):
            list[i + 1] = list[i]
        list[0] = first

    return list


