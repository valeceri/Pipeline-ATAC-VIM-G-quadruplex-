#!/usr/bin/env python3

import argparse

import pandas as pd

from Bio import SeqIO

import funz_conversione as fc

HOME = "/data/valentina.ceriani/QParse/"

FASTA = f"{HOME}fasta/"

RISULTATI = f"{HOME}risultati/"


def main(num_chr, cartella_input, regione, nome_input, min_len_loop, max_len_loop, cartella_output1, cartella_output2, cartella_output3, verso):

    # aprire file contenente la sequenza del cromosoma
    cromosoma = SeqIO.read(f"{FASTA}sequence_chr{num_chr}.fasta", "fasta")

    # aprire il file con le posizioni delle regioni di interesse
    atac_vim = pd.read_csv(f"{cartella_input}{regione}_{nome_input}.bed", sep='\t', header=None)

    # creare un DF per salvare la sequenza contenente quadruplex, il trenino di G e i loop
    treno_G = pd.DataFrame(columns=['chr', 'start', 'end', 'quadruplex'])

    treno_loop = pd.DataFrame(columns=['chr', 'start', 'end', 'pos_trenoG', 'num_loop', 'len_loop', 'loop'])

    treno_isole = pd.DataFrame(columns=['chr', 'start', 'end', 'pos_isola', 'num_isole', 'len_isole', 'isole'])

    # regioni = atac_vim.loc[atac_vim[0] == num_chr]

    # creare la lista di cromosomi
    lista_cromosomi = list(range(1, 23))
    x = ['X', 'Y', 'M']
    for el in x:
        lista_cromosomi.append(el)
    # estrarre il numero del cromosoma in ingresso da considerare per estrarre solo start e end delle regioni corrispondenti a esso
    for cr in lista_cromosomi:
        if f'chromosome {cr},' in cromosoma.description:
            cromosoma_name = f'chr{cr}'
    regioni_atac_vim = atac_vim.loc[atac_vim[0] == cromosoma_name]

    # creare file per salvare i risultati
    treno_G.to_csv(f"{cartella_output1}{verso}_{regione}_{nome_input}_treno_chr{num_chr}.bed", index=False, sep='\t')

    treno_loop.to_csv(f"{cartella_output2}{verso}_{regione}_{nome_input}_loop_chr{num_chr}.bed", index=False, sep='\t')

    treno_isole.to_csv(f"{cartella_output3}{verso}_{regione}_{nome_input}_isole_chr{num_chr}.bed", index=False, sep='\t')

    # lunghezza minima del quadruplex
    min_len = 16 + 7 * min_len_loop
    # max_len_loop = 7 (prima 12)

    # conto_treno = 0

    # per ogni sequenza compresa tra start e end controllo se ci sono quadruplex
    for indice, riga in regioni_atac_vim.iterrows():
        isola = verso + verso
        mezza_isola = verso
        pos_isola = list()
        start = riga[1]
        end = riga[2]
        pos = start
        gg = list()
        if end-start >= min_len:  # se la regione non ha la lunghezza minima va scartata
            while pos <= end:  # contare le isole
                sequenza = cromosoma[pos:end+1].seq  # estrarre la sequenza
                # print('TROVARE LE ISOLE')
                # print(sequenza)
                pos_find = sequenza.find(isola)  # con find cercare le isole, slavarle e salvare la posizione
                if pos_find != -1:
                    gg.append(isola)
                    pos = pos + pos_find
                    pos_isola.append(pos)
                    salto = 2
                    pos_find = pos_find + salto
                    altreG = True
                    contoG = 0
                    while altreG:  # se l'isola è costituita da più di 2 G l'aggiungo all'isola
                        if pos_find < len(sequenza):
                            if sequenza[pos_find] == mezza_isola:
                                pos_find += 1
                                contoG += 1
                                salto += 1
                                gg[-1] = gg[-1] + mezza_isola
                            else:
                                altreG = False
                        else:
                            altreG = False
                    pos = pos + salto
                else:
                    pos = end+1
            # print('POSIZIONE ISOLE')
            # print(pos_isola)
            # controllo loop
            if len(pos_isola) >= 8:  # se non ci sono almeno 8 isole non puà essere un doppio G4
                seq_loop = cromosoma[pos_isola[0]:pos_isola[-1]].seq
                seq_loop = str(seq_loop)
                # print('SEQUENZA PRIMA DELLA CONVERSIONE')
                # print(seq_loop)
                if isola == 'CC':  # fare l'inverso complementare
                    isola = 'GG'
                    mezza_isola = 'G'
                    seq_loop = fc.convert(seq_loop)
                    # print(gg)
                    gg = fc.convert(gg)
                    # print(gg)
                # print('SEQUENZA DOPO LA CONVERSIONE')
                # print(seq_loop)
                lista_loop = seq_loop.split(isola)
                # print(lista_loop)
                # lista_loop = [loop if loop[0] != mezza_isola else loop.replace(mezza_isola, '') for loop in lista_loop]
                # eliminare G a inizio e fine loop che vanno incluse nell'isola
                lista_loop = [loop[1:] if loop.startswith(mezza_isola) else loop for loop in lista_loop]
                lista_loop = [loop[:-1] if loop.endswith(mezza_isola) else loop for loop in lista_loop]
                lista_loop = [loop for loop in lista_loop if loop != '']
                # print(lista_loop)
                lista_loop = ['-' if len(loop) > max_len_loop else loop for loop in lista_loop]  # loop maggiore di max_len_loop non può fare quadruplex
                size = len(lista_loop) - 1
                lista_indice_loop = list(range(0, len(lista_loop)))
                sep = [idx for idx, val in enumerate(lista_loop) if val == '-']
                if len(sep) == 0:  # se non ci sono loop > max_len_loop salvo tutto
                    lista_treni = [lista_loop]
                    lista_indici_treni = [lista_indice_loop]
                else:  # sostituire il loop con - e creare liste di treni compresi tra -
                    lista_indice_loop = ['-' if i in sep else i for i in lista_indice_loop]
                    lista_treni = [lista_loop[i: j] for i, j in zip([0] + sep, sep + ([size] if sep[-1] != size else []))]
                    lista_indici_treni = [lista_indice_loop[i: j] for i, j in zip([0] + sep, sep + ([size] if sep[-1] != size else []))]
                for el in range(0, len(lista_treni)):
                    treno = lista_treni[el]
                    indice_treno = lista_indici_treni[el]
                    treno = [loop for loop in treno if loop != '-']
                    indice_treno = [it for it in indice_treno if it != '-']
                    if len(treno) >= 7:  # vedere se c'è una sequenza di almeno sette loop che abbiano lunghezza adeguata
                        # conto_treno += 1
                        quadruplex = str()
                        trenino = str()
                        isole_GG = str()
                        num_loop = len(treno)
                        len_loop = str()
                        GG = gg[indice_treno[0]:indice_treno[-1]+2]
                        num_isole = len(GG)
                        len_isole = str()
                        new_pos_isola = pos_isola[indice_treno[0]:indice_treno[-1]+2]
                        for num in range(0, len(treno)):  # mettere insieme le seq di loop, isole e treno
                            len_loop = len_loop + str(len(treno[num])) + '-'
                            trenino = trenino + treno[num] + '-'
                            len_isole = len_isole + str(len(GG[num])) + '-'
                            isole_GG = isole_GG + GG[num] + '-'
                            quadruplex = quadruplex + GG[num] + treno[num]
                        quadruplex = quadruplex + GG[-1]
                        trenino = trenino[0:-1]
                        isole_GG = isole_GG + GG[-1]
                        len_isole = len_isole + str(len(GG[-1]))
                        len_loop = len_loop[0:-1]
                        # if isola == 'CC':
                        #     quadruplex = fc.convert(quadruplex)
                        #     trenino = fc.convert(trenino)
                        #     len_loop = fc.convert(len_loop)
                        #     isole_GG = fc.convert(isole_GG)
                        #     len_isole = fc.convert(len_isole)
                        # print(quadruplex)
                        # print(trenino)
                        # print(isole_GG)
                        # salvataggio risultati
                        df_riga = pd.DataFrame({'chr': [cromosoma_name], 'start': [new_pos_isola[0]], 'end': [new_pos_isola[-1]], 'quadruplex': [quadruplex]})
                        df_loop = pd.DataFrame({'chr': [cromosoma_name], 'start': [new_pos_isola[0]], 'end': [new_pos_isola[-1]], 'pos_loop': [indice], 'num_loop': [num_loop], 'len_loop': [len_loop], 'loop': [trenino]})
                        df_isola = pd.DataFrame({'chr': [cromosoma_name], 'start': [new_pos_isola[0]], 'end': [new_pos_isola[-1]], 'pos_isola': [indice], 'num_isole': [num_isole], 'len_isole': [len_isole], 'isole': [isole_GG]})
                        # df_loop = pd.DataFrame({'pos_loop': [conto_treno], 'num_loop': [num_loop], 'len_loop': [len_loop], 'loop': [trenino]})
                        # df_isola = pd.DataFrame({'pos_isola': [conto_treno], 'num_isole': [num_isole], 'len_isole': [len_isole], 'isole': [isole_GG]})
                        df_riga.to_csv(f"{cartella_output1}{verso}_{regione}_{nome_input}_treno_chr{num_chr}.bed", mode='a', header=False, sep='\t', index=False)
                        df_loop.to_csv(f"{cartella_output2}{verso}_{regione}_{nome_input}_loop_chr{num_chr}.bed", mode='a', header=False, sep='\t', index=False)
                        df_isola.to_csv(f"{cartella_output3}{verso}_{regione}_{nome_input}_isole_chr{num_chr}.bed", mode='a', header=False, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cerca quadruplex nel cromosoma')
    parser.add_argument('num_chr', type=str)
    parser.add_argument('cartella_input', type=str)
    parser.add_argument('regione', type=str)
    parser.add_argument('nome_input', type=str)
    parser.add_argument('min_len_loop', type=int)
    parser.add_argument('max_len_loop', type=int)
    parser.add_argument('cartella_output1', type=str)
    parser.add_argument('cartella_output2', type=str)
    parser.add_argument('cartella_output3', type=str)
    parser.add_argument('verso', type=str)

    args = parser.parse_args()

    main(args.num_chr, args.cartella_input, args.regione, args.nome_input, args.min_len_loop, args.max_len_loop, args.cartella_output1, args.cartella_output2, args.cartella_output3, args.verso)

