# Needleman-Wunsch Bioinformatyka strukturalna sprawozdanie drugie.

Program zrobiony w języku python.
Instalacja:
1. Sklonowanie repozytorium
2. Instalacja bibliotek
   `pip install -r requirements.txt`
# Złożoność algorytmu
Złożoność czasowa:
Wypełnienie macierzy - O(m * n),
Traceback z algorytmem DFS - O(3^(m+n)),
gdzie:
- n - długość sekwencji pierwszej
- m - długość sekwencji drugiej

W sumie: O(n * m + 3^(n+m) * (n + m))

Złożoność pamięciowa:
Wypełnienie macierzy - O(m * n),
Traceback z algorytmem DFS - O(n+m)
gdzie:
- n - długość sekwencji pierwszej
- m - długość sekwencji drugiej

W sumie: O(n * m + k * (n + m) + (n + m))

# Przykład działania programu

Program jest możliwy do uruchomienia w trzech opcjach:
```
python needleman_wunsch.py [--fasta] file1 file2
python needleman_wunsch.py [--fasta] file1
python needleman_wunsch.py file1
python needleman_wunsch.py file1 file2
```
Argumenty pliku są obowiązkowe. Jest możliwość robienia dopasowania na dwóch plikach z wybraniem konkretnych sekwencji z pliku.

Obecne stałe w kodzie:
- match: 1
- mismatch: -1
- gap: -2

Program wykrył 3 dopasowania globalne które mają taki sam bardzo duży score **-3**. Wynik ten jest rezultatem programu - `python .\needleman_wunsch.py test.fasta`. czyli przykladu w ktorym w pliku sa tylko dwie sekwencje i automatycznie sa dopasowywane do siebie.

Oto dwie sekwencje:
```
>8V8S_1|Chains A, B, C, D|Aquaporin-4|Homo sapiens (9606)  
GATTACA  
>8V91_1|Chains A, B, C, D|Aquaporin-4|Homo sapiens (9606)  
GTCGACGCA
```
![](https://i.imgur.com/40MTQSh.png)

```
Aligning sequences of length 7 and 9...
Found 2 alignment(s) with optimal score
Optimal score: -3

Showing 2 alignment(s):

Alignment 1:
GATTA--CA
GTCGACGCA

Alignment 2:
GATTAC--A
GTCGACGCA
```
Widać na przykładzie, że alignment się zgadza.
