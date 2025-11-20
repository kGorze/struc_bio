Narzędzia Bioinformatyczne: Analiza Białek i RNA
Ten zbiór zawiera cztery skrypty w języku Python oparte na bibliotece BioPython. Służą one do analizy strukturalnej białek (mapy kontaktów, wykresy Ramachandrana) oraz do konwersji modeli struktur RNA (pomiędzy modelem pełnoatomowym a gruboziarnistym).
Wymagania
Do poprawnego działania skryptów wymagane są następujące biblioteki:
Python 3.x
Biopython
NumPy
Matplotlib
Możesz je zainstalować komendą:
pip install biopython numpy matplotlib

1. Analiza struktury Białek

contact_map.py - Mapa Kontaktów
Skrypt generuje i wyświetla mapę kontaktów dla białka 4YWO.
Działanie: Pobiera plik PDB, ekstrahuje atomy węgla alfa i oblicza macierz odległości między nimi. Jeśli odległość między atomami jest mniejsza lub równa 8 Ångströmów, punkt jest zaznaczany na mapie.
Uruchomienie: python contact_map.py

Ramachandran.py - Wykres Ramachandrana
Skrypt oblicza kąty torsyjne szkieletu białkowego i tworzy wykres Ramachandrana.
Działanie: Dla każej reszty aminokwasowej oblicza kąty Phi oraz Psi na podstawie geometrii atomów N-CA-C oraz wyświetla wykres punktowy.
Uruchomienie python Ramachandran.py

2.Przetwarzanie struktur RNA
all_to_cg.py - Konwersja do modelu gruboziarnistego
Redukuje pełną strukturę RNA do wybranych atomów reprezentatywnych (szkielet + atomy pierścienia)
    Wybrane atomy:
        Szkielet: P,C4'
        Puryny(A,G): N9,C2,C6
        Piramidyny: (C,N): N1, C2, C4
Uruchomienie: python all_to_cg.py <plik_wejściowy.pdb> <plik_wyjściowy.pdb>

cg_to_all.py - Rekonstrukcja do modelu pełnoatomowego
Odtwarza pełną strukturę atomową na podstawie modelu gruboziarnistego za pomocą szablonów.
Wymagania: Skrypt wymaga istnienia folderu templates/ w którym znajdować się będą pliki PDB dla poszczególnych nukleotydów: A.pdb, G.pdb, C.pdb, U.pdb
Działanie: Dla każdego nukleotydu w pliku wejściowym CG pobiera odpowiedni szablon, oblicza macierz rotacji i translacji (superimpozycja) na poodstawie wspólnych atomów i zapisuje zrekonstruowaną resztę.
Uruchomienie: python cg_to_all.py <plik_cg.pdb> <plik_zrekonstruowany.pdb>