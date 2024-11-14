#!/usr/bin/env python3

# QUADRUPLEX REGOLARI NEI TRENINI DI G TROVATI E CORRISPONDENTI LOOP

import argparse

import pandas as pd

HOME = "/data/valentina.ceriani/QParse/"

FASTA = f"{HOME}fasta/"


def main(cartella_input_isola, cartella_input_loop, cartella_output, verso, regione, nome_input, num_chr):
# def main(cartella_input_isola, cartella_input_loop, cartella_output, verso, regione, nome_input, num_chr, cartella_ouput_count):

    quadruplex = pd.read_csv(f"{cartella_input_isola}{verso}_{regione}_{nome_input}_isole_chr{num_chr}.bed", sep='\t')

    loop = pd.read_csv(f"{cartella_input_loop}{verso}_{regione}_{nome_input}_loop_chr{num_chr}.bed", sep='\t')

    df_out1 = pd.DataFrame(columns=['chr', 'pos_isola', 'len_isole', 'sequenza_isole', 'len_loop', 'sequenza_loop'])

    df_out1.to_csv(f"{cartella_output}{verso}_{regione}_{nome_input}_seq_corrette_chr{num_chr}.bed", sep='\t', index=False)

    # CONTARE QUANTE TRIPLETTE SENZA RIPETIZIONI
    # count_3 = 0
    # if loop.empty:
    #     df_count = pd.DataFrame([[count_3]])  # non salva l'header
    #     df_count.to_csv(f"{cartella_ouput_count}{verso}_{regione}_{nome_input}_totale_triplette_chr{num_chr}.bed", sep='\t', header=False, index=False)

    for indice, riga in quadruplex.iterrows():

        # dai file si isole e loop ottenuti con trenino estrarre le stringhe di loop e len loop, isole e len isole
        df_loop = loop.loc[indice, 'loop']
        df_loop = df_loop.split('-')

        len_loop = loop.loc[indice, 'len_loop']
        lista_loop = len_loop.split('-')

        pos_isola = riga['pos_isola']

        isole = riga['isole']
        isole = isole.split('-')

        len_isole = riga['len_isole']
        lista_isole = len_isole.split('-')

        start = 0
        end = len(lista_isole)-1
        new_end = start + 8

        # CONTARE QUANTE TRIPLETTE SENZA RIPETIZIONI
        # delta_end = new_end

        # senza ripetizioni sliding window
        lista_loop_df = list()

        # con sliding window selezionare 8 isole alla volta e controllare che per entrambi i G4 le isole centrali siano
        # uguali e di 2 o 3 G, mentre quelle laterali siano di 2 (se quelle centrali sono da 2), 3 o 4. Salvare loop,
        # isole e lunghezze corrispondenti al quadruplex doppio
        while new_end <= end:

            new_lista = lista_isole[start:new_end]
            new_isole = isole[start:new_end]
            new_loop = lista_loop[start:new_end]
            new_df_loop = df_loop[start:new_end]

            isole_centrali = [lista_isole[start+1], lista_isole[start+2], lista_isole[start+5], lista_isole[start+6]]
            isole_esterne = [lista_isole[start], lista_isole[start+3], lista_isole[start+4], lista_isole[start+7]]

            # condizione per isole centrali: uguali nel primo e nel secondo g4 e lunghezza inferiore a 4
            if isole_centrali[0] == isole_centrali[1] and isole_centrali[2] == isole_centrali[3] and all(int(y) < 4 for y in isole_centrali):

                # condizione per isole laterali: devono essere di lunghezza maggiore o uguale alle isole centrali e inferiore a 5
                if all(int(x) <= 4 for x in isole_esterne) and all(int(a) <= int(b) for a, b in zip(isole_centrali, isole_esterne)):

                    len_isole_df = str()
                    isole_df = str()
                    len_loop_df = str()
                    loop_df = str()

                    for el in range(0, 7):
                        len_isole_df = len_isole_df + new_lista[el] + '-'
                        isole_df = isole_df + new_isole[el] + '-'
                        len_loop_df = len_loop_df + new_loop[el] + '-'
                        loop_df = loop_df + new_df_loop[el] + '-'

                    len_isole_df = len_isole_df + new_lista[-1]
                    isole_df = isole_df + new_isole[-1]

                    len_loop_df = len_loop_df.rstrip('-')
                    loop_df = loop_df.rstrip('-')


                    


                    # senza ripetizioni sliding window
                    if loop_df not in lista_loop_df:
                        lista_loop_df.append(loop_df)
                        df_q = pd.DataFrame({'chr': [num_chr], 'pos_isola': [pos_isola], 'len_isole': [len_isole_df], 'sequenza_isole': [isole_df], 'len_loop': [len_loop_df], 'sequenza_loop': [loop_df]})
                        df_q.to_csv(f"{cartella_output}{verso}_{regione}_{nome_input}_seq_corrette_chr{num_chr}.bed", mode='a', header=False, sep='\t', index=False)
                        start += 1
                        new_end += 1
                    else:
                        start += 1
                        new_end += 1

                    # originale
                    # df_q = pd.DataFrame({'chr': [num_chr], 'pos_isola': [pos_isola], 'len_isole': [len_isole_df], 'sequenza_isole': [isole_df], 'len_loop': [len_loop_df], 'sequenza_loop': [loop_df]}) # lista_loop_df
                    # df_q.to_csv(f"{cartella_output}{verso}_{regione}_{nome_input}_seq_corrette_chr{num_chr}.bed", mode='a', header=False, sep='\t', index=False)
                    # start += 1
                    # new_end += 1

                    # CONTARE QUANTE TRIPLETTE SENZA RIPETIZIONI
                    # if delta_end - start == 0:
                    #     delta_end = new_end
                    #     lista_len = len_loop_df.split('-')
                    #     lista_int = [int(x) for x in lista_len]
                    #     conteggio = lista_int.count(3)
                    #     count_3 = count_3 + conteggio

                else:
                    idx = [index for index, value in enumerate(new_lista) if int(value) > 4]
                    if len(idx) == 0:
                        idx = [0]
                    start = start + idx[-1] + 1
                    new_end = start + 8

            else:
                start += 1
                new_end += 1

            # CONTARE QUANTE TRIPLETTE SENZA RIPETIZIONI
            # df_count = pd.DataFrame([[count_3]])  # non salva l'header
            # df_count.to_csv(f"{cartella_ouput_count}{verso}_{regione}_{nome_input}_totale_triplette_chr{num_chr}.bed", sep='\t', header=False, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Lunghezze dei loop in colonna')
    parser.add_argument('cartella_input_isola', type=str)
    parser.add_argument('cartella_input_loop', type=str)
    parser.add_argument('cartella_output', type=str)
    parser.add_argument('verso', type=str)
    parser.add_argument('regione', type=str)
    parser.add_argument('nome_input', type=str)
    parser.add_argument('num_chr', type=str)
    # parser.add_argument('cartella_ouput_count', type=str)

    args = parser.parse_args()

    main(args.cartella_input_isola, args.cartella_input_loop, args.cartella_output, args.verso, args.regione, args.nome_input, args.num_chr)
    # main(args.cartella_input_isola, args.cartella_input_loop, args.cartella_output, args.verso, args.regione, args.nome_input, args.num_chr, args.cartella_ouput_count)

