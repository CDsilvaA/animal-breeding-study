#!/bin/sh

# Arquivos obrigatórios (exemplos)
# data.txt
# 12000 0 2 1 0 5 0
# 12001 0 2 1 0 5 0
# 12002 0 2 1 0 5 0

# Lista das colunas para serem selecionadas
# list.txt
# 2
# 3
# 5

# Arquivo de texto
# ID.txt
# ID

# Rodar comandos na ordem:

# informa o total de colunas dos dados
awk '{print NF}' data.txt | uniq

# cria a sequência das colunas (vetor coluna)
# seq 1 (total_colunas - 1) > seq
seq 1 6 > seq

# transpõe a sequência (vetor linha)
awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' seq | tr ':' ' ' >> temp1

# junta o arquivo ID com o arquivo de sequencia transposto
paste -d '' ID temp1 >> temp2

# adiciona o cabeçalho no arquivo de dados
printf '0r temp2\nwq\n' | ed -s data.txt

# seleciona as colunas a partir de uma lista externa
cols=($(sed '1!d;s/ /\n/g' data.txt | grep -nf list.txt | sed 's/:.*$//'))
cut -d ' ' -f 1$(printf ",%s" "${cols[@]}") data.txt >> temp.txt
cat temp.txt 
sed 1d temp.txt > output.txt
rm temp.txt temp1 temp2 seq.txt

# adiciona cabeçalho (se precisar)
sed -e '1i\HEADER_HERE' data > data.txt

# Fonte: https://stackoverflow.com/questions/11098189/extract-columns-from-a-file-based-on-header-selected-from-another-file
# Fonte: https://unix.stackexchange.com/questions/25592/creating-a-sequence-of-numbers-one-per-line-in-a-file
# Fonte: https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash

