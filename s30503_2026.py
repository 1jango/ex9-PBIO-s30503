import random

# Numer albumu: s30503
# Data: 2026-05-12
# Opis: Generator sekwencji DNA - wersja podstawowa.

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Pobiera liczbę od użytkownika z walidacją."""
    while True:
        try:
            val = int(input(prompt))
            if min_val <= val <= max_val: return val
            print(f"Error: range [{min_val}, {max_val}]")
        except ValueError:
            print("Error: enter an integer.")

def generate_sequence(length: int) -> str:
    """Generuje losowe DNA."""
    return "".join(random.choice("ACGT") for _ in range(length))

def main():
    length = validate_positive_int("Enter sequence length: ")
    raw_dna = generate_sequence(length)
    print(f"Generated {len(raw_dna)} nucleotides.")

if __name__ == "__main__":
    main()