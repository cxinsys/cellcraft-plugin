import sys
import pandas as pd
import plotly.express as px
import json
import traceback

def log_error(message, exception=None):
    """Print error message with optional exception details."""
    print(f"[ERROR] {message}", file=sys.stderr)
    if exception:
        print(f"[ERROR] Exception: {type(exception).__name__}: {str(exception)}", file=sys.stderr)
        print(f"[ERROR] Traceback:\n{traceback.format_exc()}", file=sys.stderr)

def log_info(message):
    """Print informational message."""
    print(f"[INFO] {message}")

def main():
    # Validate command line arguments
    if len(sys.argv) < 4:
        log_error(f"Insufficient arguments. Expected 3 arguments (input_file, output_file, top_n), got {len(sys.argv) - 1}")
        log_error("Usage: python Barplot.py <input_file> <output_file> <top_n>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    # Validate top_n argument
    try:
        top_n = int(sys.argv[3])
        if top_n <= 0:
            log_error(f"top_n must be a positive integer, got {top_n}")
            sys.exit(1)
    except ValueError as e:
        log_error(f"Invalid top_n value: '{sys.argv[3]}'. Must be an integer.", e)
        sys.exit(1)

    log_info(f"Processing input file: {input_file_path}")
    log_info(f"Output file: {output_file_path}")
    log_info(f"Top N genes: {top_n}")

    # Load data and check format
    try:
        df = pd.read_csv(input_file_path, sep=None, engine='python', header=None)
        log_info(f"Successfully loaded data with {len(df)} rows and {len(df.columns)} columns")
    except FileNotFoundError as e:
        log_error(f"Input file not found: {input_file_path}", e)
        sys.exit(1)
    except pd.errors.EmptyDataError as e:
        log_error(f"Input file is empty: {input_file_path}", e)
        sys.exit(1)
    except Exception as e:
        log_error(f"Failed to read input file: {input_file_path}", e)
        sys.exit(1)

    # Validate data is not empty
    if df.empty:
        log_error("Input file contains no data rows")
        sys.exit(1)

    # Process data based on column count
    try:
        if len(df.columns) == 3:  # fdrTrim.sif format (3 columns)
            df.columns = ['Source', 'Weight', 'Target']
            # Calculate Target count per Source
            source_freq = df['Source'].value_counts().reset_index()
            source_freq.columns = ['Genes', 'Outbound']
            log_info("Detected 3-column format (Source, Weight, Target)")
        elif len(df.columns) == 2:  # numlinksOutdegree.txt format (2 columns)
            df.columns = ['Genes', 'Outbound']
            source_freq = df.copy()
            log_info("Detected 2-column format (Genes, Outbound)")
        else:
            log_error(f"Input file must have 2 or 3 columns, got {len(df.columns)} columns")
            sys.exit(1)
    except Exception as e:
        log_error("Failed to process data columns", e)
        sys.exit(1)

    # Select top N genes
    actual_top_n = min(top_n, len(source_freq))
    if actual_top_n < top_n:
        log_info(f"Requested top {top_n} genes, but only {actual_top_n} available")
    source_freq = source_freq.head(actual_top_n)

    # Sort in descending order
    source_freq = source_freq.sort_values(by='Outbound', ascending=False)

    # Custom color palette
    color_palette = px.colors.qualitative.Set3

    # Create visualization with Plotly
    try:
        fig = px.bar(
            source_freq,
            x='Outbound',
            y='Genes',
            color='Genes',
            orientation='h',
            title=f"Top {actual_top_n} Genes Outbound Barplot",
            color_discrete_sequence=color_palette
        )

        fig.update_layout(yaxis={'categoryorder':'total ascending'})
        log_info("Successfully created barplot visualization")
    except Exception as e:
        log_error("Failed to create Plotly visualization", e)
        sys.exit(1)

    # Save as JSON
    try:
        fig_json = fig.to_json()
        with open(output_file_path, "w") as json_file:
            json_file.write(fig_json)
        log_info(f"Successfully saved output to: {output_file_path}")
    except PermissionError as e:
        log_error(f"Permission denied when writing to: {output_file_path}", e)
        sys.exit(1)
    except Exception as e:
        log_error(f"Failed to save output file: {output_file_path}", e)
        sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_error("Unexpected error occurred", e)
        sys.exit(1)
