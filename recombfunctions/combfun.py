from random import randint


# Chooses a random individual from a list of genetic combinations
def rand_ind(comb_list):
    random_individual = comb_list[randint(0, len(comb_list) - 1)]
    return random_individual


# Chooses randomly one of the two EF1A variants of an individual
def choose_allele(ef1a):
    allele = ef1a[randint(0, 1)]
    return allele


# Generates a new individual (combination of genetic variants) from genetic variants of two parents
def combine(ind1, ind2):
    randssu = (ind1[0], ind2[0])[randint(0, 1)]
    randcoi = (ind1[2], ind2[2])[randint(0, 1)]
    randef = [choose_allele(ind1[1]), choose_allele(ind2[1])]
    return [randssu, randef, randcoi]


# Generates a new individual with unequal weights for SSU and COI variants
def combine_unequal(ind1, ind2, weights_ssu, weights_coi):
    # choose SSU according to its weight
    if weights_ssu[ind1[0]] > weights_ssu[ind2[0]]:
        randssu = ind1[0]
    elif weights_ssu[ind2[0]] > weights_ssu[ind1[0]]:
        randssu = ind2[0]
    else:
        randssu = (ind1[0], ind2[0])[randint(0, 1)]
    # choose CI according to its weight
    if weights_coi[ind1[2]] > weights_coi[ind2[2]]:
        randcoi = ind1[2]
    elif weights_coi[ind2[2]] > weights_coi[ind1[2]]:
        randcoi = ind2[2]
    else:
        randcoi = (ind1[2], ind2[2])[randint(0, 1)]
    # choose EF1A (Mendelian)
    randef = [choose_allele(ind1[1]), choose_allele(ind2[1])]
    return [randssu, randef, randcoi]


def write_log(file, string):
    with open(file, "a") as f:
        print(string)
        f.write(string + "\n")
