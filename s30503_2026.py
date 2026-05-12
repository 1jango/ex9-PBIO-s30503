import random


# s30503, 2026-05-12
# Etap 2: Dodanie obsługi formatu FASTA, obliczanie zawartości GC i wstawianie imienia.

def validate_positive_int(prompt, min_val=1, max_val=100_000):
    """Pobiera liczbę od użytkownika i sprawdza czy mieści się w zakresie."""
    while True:
        try:
            line = input(prompt)
            value = int(line)
            if min_val <= value <= max_val:
                return value
            print(f"Błąd: Wartość musi być w przedziale {min_val}-{max_val}.")
        except ValueError:
            print("Błąd: To nie jest poprawna liczba całkowita.")


def generate_sequence(length: int) -> str:
    """Generuje losowy ciąg nukleotydów DNA."""
    zasady = "ACGT"
    return "".join(random.choice(zasady) for _ in range(length))


def calculate_stats(sequence: str) -> dict:
    """Oblicza procentowy skład nukleotydów. Ignoruje małe litery (imię)."""
    # Wyciągamy tylko wielkie litery, żeby statystyki były rzetelne
    dna_only = "".join([c for c in sequence if c.isupper()])
    n = len(dna_only)

    if n == 0:
        return {b: 0.0 for b in "ACGT"} | {"GC": 0.0}

    stats = {}
    for base in "ACGT":
        count = dna_only.count(base)
        stats[base] = (count / n) * 100

    # Standardowe obliczenie zawartości GC
    stats["GC"] = stats["G"] + stats["C"]
    return stats


def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię (małymi literami) w losowe miejsce sekwencji."""
    index = random.randint(0, len(sequence))
    return sequence[:index] + name.lower() + sequence[index:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Tworzy rekord w formacie FASTA z zawijaniem wierszy."""
    header = f">{seq_id} {description}".strip()
    # Dzielenie sekwencji na linie o stałej szerokości
    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]
    return header + "\n" + "\n".join(lines) + "\n"


def main():
    # Pobieranie danych wejściowych
    length = validate_positive_int("Podaj długość sekwencji: ")
    seq_id = input("Podaj ID sekwencji: ").strip().replace(" ", "_")
    desc = input("Podaj opis sekwencji: ")
    user_name = input("Podaj swoje imię: ")

    # Logika generowania
    raw_dna = generate_sequence(length)
    stats = calculate_stats(raw_dna)

    # Przygotowanie finalnej sekwencji z wstawionym imieniem
    dna_with_name = insert_name(raw_dna, user_name)

    # Zapis do pliku
    filename = f"{seq_id}.fasta"
    with open(filename, "w") as f:
        f.write(format_fasta(seq_id, desc, dna_with_name))

    # Wyświetlanie wyników
    print(f"\nSekwencja została zapisana w pliku: {filename}")
    print(f"Statystyki składu (n={length}):")
    for base in "ACGT":
        print(f"  {base}: {stats[base]:.2f}%")
    print(f"GC-content: {stats['GC']:.2f}%")


if __name__ == "__main__":
    main()