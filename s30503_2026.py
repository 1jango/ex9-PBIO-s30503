import random


# s30503, 2026-05-12
# Projekt: Zaawansowany generator i analizator sekwencji DNA w formacie FASTA.

def sprawdz_liczbe(pytanie, dolna=1, gorna=100_000):
    """Pobiera liczbe od uzytkownika i dba o poprawnosc zakresu."""
    while True:
        try:
            n = int(input(pytanie))
            if dolna <= n <= gorna:
                return n
            print(f"Blad: Liczba musi byc w przedziale {dolna}-{gorna}")
        except ValueError:
            print("Blad: To nie jest liczba calkowita.")


def stworz_dna(dlugosc):
    """Tworzy losowy ciag nukleotydow (tylko wielkie litery)."""
    return "".join(random.choice("ACGT") for _ in range(dlugosc))


def analizuj_sklad(sekwencja):
    """Oblicza procentowa zawartosc nukleotydow i GC-content."""
    # Filtrujemy, aby liczyc tylko DNA (wielkie litery), pomijajac imie autora
    czyste_dna = "".join([znak for znak in sekwencja if znak.isupper()])
    ile = len(czyste_dna)

    if ile == 0:
        return {b: 0.0 for b in "ACGT"} | {"GC": 0.0}

    wyniki = {b: (czyste_dna.count(b) / ile) * 100 for b in "ACGT"}
    wyniki["GC"] = wyniki["G"] + wyniki["C"]
    return wyniki


def wstaw_imie(dna, imie):
    """Wstawia imie autora (malymi literami) w losowe miejsce sekwencji."""
    punkt = random.randint(0, len(dna))
    return dna[:punkt] + imie.lower() + dna[punkt:]


def zapisz_fasta(id_s, opis, dna, szerokosc=80):
    """Zwraca sformatowany rekord FASTA z zawijaniem wierszy."""
    naglowek = f">{id_s} {opis}".strip()
    # Rozbijanie sekwencji na linie o szerokosci 80 znakow
    linie = [dna[i:i + szerokosc] for i in range(0, len(dna), szerokosc)]
    return naglowek + "\n" + "\n".join(linie) + "\n"


# --- DODATKOWE FUNKCJONALNOSCI BIOINFORMATYCZNE ---

def transkrypcja(dna):
    """Proces zamiany DNA na mRNA (T zamieniane na U)."""
    return dna.replace("T", "U")


def znajdz_motyw(dna, motyw):
    """Szuka pozycji podanego motywu w sekwencji (indeksowanie od 1)."""
    dna_czyste = "".join([c for c in dna if c.isupper()])
    pozycje = []
    idx = dna_czyste.find(motyw.upper())
    while idx != -1:
        pozycje.append(idx + 1)  # Biologiczne indeksowanie zaczyna sie od 1
        idx = dna_czyste.find(motyw.upper(), idx + 1)
    return pozycje


def odwroc_i_uzupelnij(dna):
    """Generuje sekwencje odwrotnie komplementarna."""
    pary = str.maketrans("ACGTacgt", "TGCAtgca")
    komplementarna = dna.translate(pary)
    return komplementarna[::-1]


def tlumacz_na_bialko(dna):
    """Tlumaczy kodony DNA na sekwencje aminokwasowa."""
    tabela = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_', 'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    czyste = "".join([c for c in dna if c.isupper()])
    bialko = ""
    # Przetwarzamy trojkami (kodonami)
    for i in range(0, len(czyste) - 2, 3):
        kodon = czyste[i:i + 3]
        bialko += tabela.get(kodon, "?")
    return bialko


def main():
    # 1. Interakcja z uzytkownikiem
    rozmiar = sprawdz_liczbe("Wprowadz dlugosc sekwencji (1-100000): ")
    identyfikator = input("Wprowadz ID (bez spacji): ").strip().replace(" ", "_")
    opis_sek = input("Wprowadz opis (moze byc pusty): ")
    imie_autora = input("Wprowadz swoje imie: ")

    # 2. Generowanie danych bazowych
    kod_dna = stworz_dna(rozmiar)
    staty = analizuj_sklad(kod_dna)
    dna_z_imieniem = wstaw_imie(kod_dna, imie_autora)

    # 3. Zapis do pliku wielo-FASTA (Multi-FASTA)
    sciezka = f"{identyfikator}.fasta"
    with open(sciezka, "w") as plik:
        # Glowny rekord DNA
        plik.write(zapisz_fasta(identyfikator, opis_sek, dna_z_imieniem))

        # Dodatkowy rekord mRNA
        rna = transkrypcja(dna_z_imieniem)
        plik.write(zapisz_fasta(f"{identyfikator}_mRNA", "Sekwencja transkrybowana", rna))

        # Dodatkowy rekord Reverse Complement
        rev = odwroc_i_uzupelnij(dna_z_imieniem)
        plik.write(zapisz_fasta(f"{identyfikator}_Reverse", "Nic odwrotnie komplementarna", rev))

        # Dodatkowy rekord bialka
        aa = tlumacz_na_bialko(kod_dna)
        plik.write(zapisz_fasta(f"{identyfikator}_Protein", "Sekwencja aminokwasowa", aa))

    # 4. Wyswietlanie statystyk i interakcja końcowa
    print(f"\nSukces! Wszystkie rekordy zapisano w pliku: {sciezka}")
    print("-" * 30)
    print(f"Statystyki dla n={rozmiar}:")
    for nukleotyd, wartosc in staty.items():
        if nukleotyd != "GC":
            print(f"  {nukleotyd}: {wartosc:.2f}%")
    print(f"Zawartosc GC: {staty['GC']:.2f}%")

    # Szukanie motywu na prosbe uzytkownika
    szukany = input("\nChcesz wyszukac motyw (np. ATG)? Podaj go lub wcisnij Enter: ")
    if szukany:
        znalezione = znajdz_motyw(kod_dna, szukany)
        print(f"Motyw '{szukany}' znaleziono na pozycjach: {znalezione}")


if __name__ == "__main__":
    main()