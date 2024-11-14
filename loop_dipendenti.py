#!/usr/bin/env python3

import argparse

import pandas as pd


HOME = "/data/valentina.ceriani/QParse/"

RISULTATI = f"{HOME}risultati/"


# G4 per G4 controllo lunghezze e sequenze verificando la posizione in cui si trovano e la frequenza totale


def main(cartella_input_loop, cartella_input_seq, regione, nome_input, cartella_output_loop, cartella_output_seq):

    # LUNGHEZZE

    df_in_len = pd.read_csv(f"{cartella_input_loop}{regione}_{nome_input}_loop_lunghezze.bed", sep='\t', dtype={0: str})

    df_out_len = pd.DataFrame(columns=['names', 'values_loop1', 'values_loop1_norm', 'values_loop2', 'values_loop2_norm', 'values_loop3', 'values_loop3_norm', 'values_tot', 'values_tot_norm'])

    df_out_len.to_csv(f"{cartella_output_loop}{regione}_{nome_input}_lunghezze_loop_dipendenti.bed", sep='\t', index=False)

    loop = list(df_in_len[f"-{regione}-"])

    loop = [el for el in loop if el != f"-{regione}-"]

    x = 0

    lista_loop1 = list()
    lista_loop2 = list()
    lista_loop3 = list()

    while x < len(loop):

        # print('x', x)

        loop1 = loop[x]
        loop2 = loop[x+1]
        loop3 = loop[x+2]

        # print(loop1, loop2, loop3)

        lista_loop1.append(int(loop1))
        lista_loop2.append(int(loop2))
        lista_loop3.append(int(loop3))

        x = x + 3

    names_loop = list()
    values_l1 = list()
    values_l2 = list()
    values_l3 = list()

    # presi primo, secondo e terzo loop separatamente li assegno alla lista corrispondente
    # prendo ogni valore una sola volta e poi conto quante volte è presente nelle liste di loop
    # calcolo la frequenza totale e i valori normalizzati (uguale per le sequenze)
    for i in range(1, 13):
        names_loop.append(i)
        # print('names_loop', names_loop)
        val1 = lista_loop1.count(i)
        # print(val1)
        values_l1.append(int(val1))
        # print('values_l1', values_l1)
        val2 = lista_loop2.count(i)
        values_l2.append(int(val2))
        # print('values_l2', values_l2)
        val3 = lista_loop3.count(i)
        values_l3.append(int(val3))
        # print('values_l3', values_l3)

    len_tot = [int(a) + int(b) + int(c) for a, b, c in zip(values_l1, values_l2, values_l3)]

    somma_val1 = sum(values_l1)
    somma_val2 = sum(values_l2)
    somma_val3 = sum(values_l3)
    somma_tot = sum(len_tot)

    values_l1_norm = [val1_norm / somma_val1 for val1_norm in values_l1]
    values_l2_norm = [val2_norm / somma_val2 for val2_norm in values_l2]
    values_l3_norm = [val3_norm / somma_val3 for val3_norm in values_l3]
    len_tot_norm = [len_norm / somma_tot for len_norm in len_tot]

    # print('len_tot', len_tot)

    df_loop = pd.DataFrame()

    df_loop['names'] = names_loop
    df_loop['values_loop1'] = values_l1
    df_loop['values_loop1_norm'] = values_l1_norm
    df_loop['values_loop2'] = values_l2
    df_loop['values_loop2_norm'] = values_l2_norm
    df_loop['values_loop3'] = values_l3
    df_loop['values_loop3_norm'] = values_l3_norm
    df_loop['values_tot'] = len_tot
    df_loop['values_tot_norm'] = len_tot_norm

    df_loop.to_csv(f"{cartella_output_loop}{regione}_{nome_input}_lunghezze_loop_dipendenti.bed", mode='a', header=False, sep='\t', index=False)

    # SEQUENZE

    df_in_seq = pd.read_csv(f"{cartella_input_seq}{regione}_{nome_input}_loop_basi_lunghezze.bed", sep='\t', dtype={0: str})

    df_out_seq = pd.DataFrame(columns=['names', 'values_seq1', 'values_seq1_norm', 'values_seq2', 'values_seq2_norm', 'values_seq3', 'values_seq3_norm', 'values_tot', 'values_tot_norm'])

    df_out_seq.to_csv(f"{cartella_output_seq}{regione}_{nome_input}_sequenze_loop_dipendenti.bed", sep='\t', index=False)

    seq = list(df_in_seq[f"-{regione}-"])

    seq = [el for el in seq if el != f"-{regione}-"]

    y = 0

    lista_seq = list()
    lista_seq1 = list()
    lista_seq2 = list()
    lista_seq3 = list()

    while y < len(seq):

        # print('y', y)

        seq1 = seq[y]
        seq2 = seq[y+1]
        seq3 = seq[y+2]

        lista_seq1.append(seq1)
        lista_seq2.append(seq2)
        lista_seq3.append(seq3)

        # print(seq1, seq2, seq3)

        y = y + 3

        if seq1 not in lista_seq:
            lista_seq.append(seq1)

        if seq2 not in lista_seq:
            lista_seq.append(seq2)

        if seq3 not in lista_seq:
            lista_seq.append(seq3)

    names_seq = list()
    values_s1 = list()
    values_s2 = list()
    values_s3 = list()

    for i in lista_seq:
        names_seq.append(i)
        # print('names_seq', names_seq)
        val1 = lista_seq1.count(i)
        values_s1.append(int(val1))
        # print('values_s1', values_s1)
        val2 = lista_seq2.count(i)
        values_s2.append(int(val2))
        # print('values_s2', values_s2)
        val3 = lista_seq3.count(i)
        values_s3.append(int(val3))
        # print('values_s3', values_s3)

    seq_tot = [int(a) + int(b) + int(c) for a, b, c in zip(values_s1, values_s2, values_s3)]

    # print('seq_tot', seq_tot)

    somma_seq1 = sum(values_s1)
    somma_seq2 = sum(values_s2)
    somma_seq3 = sum(values_s3)
    somma_seq_tot = sum(seq_tot)

    values_s1_norm = [seq1_norm / somma_seq1 for seq1_norm in values_s1]
    values_s2_norm = [seq2_norm / somma_seq2 for seq2_norm in values_s2]
    values_s3_norm = [seq3_norm / somma_seq3 for seq3_norm in values_s3]
    seq_tot_norm = [seq_norm / somma_seq_tot for seq_norm in seq_tot]

    df_seq = pd.DataFrame()

    df_seq['names'] = names_seq
    df_seq['values_seq1'] = values_s1
    df_seq['values_seq1_norm'] = values_s1_norm
    df_seq['values_seq2'] = values_s2
    df_seq['values_seq2_norm'] = values_s2_norm
    df_seq['values_seq3'] = values_s3
    df_seq['values_seq3_norm'] = values_s3_norm
    df_seq['values_tot'] = seq_tot
    df_seq['values_tot_norm'] = seq_tot_norm

    df_seq.to_csv(f"{cartella_output_seq}{regione}_{nome_input}_sequenze_loop_dipendenti.bed", mode='a', header=False, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Grafico distribuzione lunghezze loop')
    parser.add_argument('cartella_input_loop', type=str)
    parser.add_argument('cartella_input_seq', type=str)
    parser.add_argument('regione', type=str)
    parser.add_argument('nome_input', type=str)
    parser.add_argument('cartella_output_loop', type=str)
    parser.add_argument('cartella_output_seq', type=str)

    args = parser.parse_args()

    main(args.cartella_input_loop, args.cartella_input_seq, args.regione, args.nome_input, args.cartella_output_loop, args.cartella_output_seq)
