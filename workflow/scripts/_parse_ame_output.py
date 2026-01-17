import sys
import io

def process_data():
    """
    Reads data line-by-line from stdin, splits fields using variable whitespace,
    replaces the first field (original rank) with the line number, and prints
    the result to stdout as clean Tab-Separated Values (TSV).
    """
    # Use TextIOWrapper for robust handling of standard input
    input_stream = io.TextIOWrapper(sys.stdin.buffer, encoding='utf-8')
    row_number = 1

    # Define the output separator as a single tab
    OUTPUT_DELIMITER = '\t'

    for line in input_stream:
        # 1. Split on any sequence of whitespace (spaces or tabs).
        # This is Python's equivalent of awk's default field separation logic.
        parts = line.strip().split()

        # Skip empty lines
        if not parts:
            continue

        # 2. Build the new row:
        #    - Prepend the current row number (as a string).
        #    - Append all parts starting from the second element (parts[1:]),
        #      effectively dropping the original first column (parts[0]).
        new_row_parts = [str(row_number)] + parts[1:]

        # 3. Join the new parts with the consistent tab delimiter and print.
        print(OUTPUT_DELIMITER.join(new_row_parts))

        row_number += 1

if __name__ == '__main__':
    process_data()