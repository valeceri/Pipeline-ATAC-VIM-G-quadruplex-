#!/usr/bin/env python3

import pandas as pd


def convert(seq):
    # input: sequenza di cui fare l'inverso complementare

    if type(seq) == str:  # se la sequenza è una stringa

        seq_inv = seq[::-1]  # invertirla

        seq_compl = str()

        # per ogni base sostituire il complementare
        for el in seq_inv:
            if el == 'A':
                seq_compl = seq_compl + 'T'
            elif el == 'G':
                seq_compl = seq_compl + 'C'
            elif el == 'T':
                seq_compl = seq_compl + 'A'
            elif el == 'C':
                seq_compl = seq_compl + 'G'
            else:
                seq_compl = seq_compl + el

    else:  # se lista di isole
        seq.reverse()  # invertire
        seq_compl = list()
        for el in seq:
            seq_compl.append(len(el)*'G')  # trasformare nel complementare

    return seq_compl




    # seq_inv = seq[::-1]
    #
    # seq_compl = str()
    #
    # for el in seq_inv:
    #     if el == 'A':
    #         seq_compl = seq_compl + 'T'
    #     elif el == 'G':
    #         seq_compl = seq_compl + 'C'
    #     elif el == 'T':
    #         seq_compl = seq_compl + 'A'
    #     elif el == 'C':
    #         seq_compl = seq_compl + 'G'
    #     else:
    #         seq_compl = seq_compl + el
    # return seq_compl