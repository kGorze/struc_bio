# Narzędzia Bioinformatyczne: Analiza Białek i RNA

Programy zostały wykonane w w języku `Python`, oparte na bibliotece `BioPython`. Służą one do analizy strukturalnej białek (mapy kontaktów, wykresy Ramachandrana) oraz do konwersji modeli struktur RNA (pomiędzy modelem pełnoatomowym a gruboziarnistym).
# Wymagania

Do poprawnego działania skryptów wymagane są następujące biblioteki:
`Python 3.x`, `Biopython`, `NumPy`, `Matplotlib`

Można zainstalować je pipem za pomocą komendy: `pip install biopython numpy matplotlib`
# Zadanie 1/2 
### Analiza struktury Białek

1. **Skrypt `contact_map.py` generuje i wyświetla mapę kontaktów dla białka** `4YWO`.
Uruchomienie: `python contact_map.py`


![](https://i.imgur.com/gjb8Ed4.png)
***Fig. 1***
*Przykładowa mapa kontaktów dla białka `4YWO`.*

2. **Działanie**
Pobiera plik PDB, ekstrahuje atomy węgla alfa i oblicza macierz odległości między nimi. Jeśli odległość między atomami jest mniejsza lub równa 8 Ångströmów, punkt jest zaznaczany na mapie.

---

1. **Skrypt `Ramachandran.py` oblicza kąty torsyjne szkieletu białkowego i tworzy wykres Ramachandrana  dla białka** `4YWO`.
Uruchomienie: `python Ramachandran.py`


![](https://i.imgur.com/FsqNTEz.png)

**Fig. 2**
*Przykładowy wykres Ramachandran dla białka `4YWO`.*

2. **Działanie**
 Dla każej reszty aminokwasowej oblicza kąty torsyjne $\Phi$ oraz $\psi$ na na podstawie czterech kolejnych atomów:  $\phi$ z $C(i-1)-N(i)-C \alpha(i) -C(i)$ oraz  $\psi$ z $N(i)-C\alpha(i)-C(i)-N(i+1)$, a następnie wyświetla wykres punktowy (wykres Ramachandrana).


# Zadanie 3 

### Przetwarzanie struktur RNA

1. **Skrypt `all_to_cg.py ` redukuje pełną strukturę RNA do wybranych atomów reprezentatywnych (szkielet + atomy pierścienia) dla wybranego pliku.**
Uruchomienie: `python all_to_cg.py <plik_wejściowy.pdb> <plik_wyjściowy.pdb>`
	
Atomy które zostały wybrane do reprezentacji coarse grained:
```
Szkielet: P,C4
Puryny(A,G): N9,C2,C6
Piramidyny: (C,N): N1, C2, C4
```


![](https://i.imgur.com/K4xxbPi.png)
![](https://i.imgur.com/cw9y5aO.png)
![](https://i.imgur.com/JRGNPVM.png)


**Fig. 3/4**
*Sturuktura grubo ziarnista nałożona na pełną. Kolor czerwony oznacza kulki które sią gruboziarniste w porównaniu do szarej pełnej struktury. Jest to RNA z pliku `430D` i `430D-gc`.*

2. **Działanie**
Skrypt wykorzystuje klasę `PDB.Select` z biblioteki *BioPython* do filtrowania atomów podczas zapisu. Z oryginalnej struktury zachowywane są wyłącznie atomy zdefiniowane jako reprezentatywne dla modelu uproszczonego.

---

1. **Skrypt `cg_to_all.py ` rekonstrułuje pełną strukturę RNA za pomocą szablonów dla wybranego pliku do modelu pełnoatomowego.**
Uruchomienie: `python cg_to_all.py <plik_cg.pdb> <plik_zrekonstruowany.pdb>`

Wymagania: ***Skrypt wymaga istnienia folderu templates/ w którym znajdować się będą pliki PDB dla poszczególnych nukleotydów: A.pdb, G.pdb, C.pdb, U.pdb***

![](https://i.imgur.com/k24qQYI.png)


**Fig. 5**
*Stuktura oryginalna i nałożona z widocznym wyliczonym `RMSD`, które wskazuje na prawie idealne odwzorowanie. Jest to RNA z pliku `430D-rc` i `430D`.

2. **Działanie**
Dla każdego nukleotydu w pliku wejściowym CG pobiera odpowiedni szablon, oblicza macierz rotacji i translacji (superimpozycja) na poodstawie wspólnych atomów i zapisuje zrekonstruowaną resztę.

3. **Wniosek**:
RMSD, które nie jest równe zeru - jednak jest jemu bliskie - wynika z utraty informacji o elastyczności i subtelnych różnicach konformacynych poszczególnych nukleotydow podczas konwersji do modelu gruboziarnistego. 
Jest to naturalny koszt upraszczania modelu CG. 