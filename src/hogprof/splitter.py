import sys
import os
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def add_fastoma_path(fastoma_path):
    """Add the FastOMA utils directory to Python's sys.path."""
    utils_path = os.path.join(os.path.abspath(fastoma_path), "utils")

    if not os.path.isdir(utils_path):
        raise FileNotFoundError(
            f"FastOMA utils directory not found: {utils_path}"
        )

    if utils_path not in sys.path:
        sys.path.append(utils_path)

def main():
    parser = argparse.ArgumentParser(description="Split OrthoXML file by species.")
    parser.add_argument("--orthoxml_dir", required=True, help="Path to the directory containing FastOMA_HOGs.orthoxml")
    parser.add_argument(
        "--fastoma-path",
        required=True,
        help="Path to the FastOMA installation directory",
    )
    args = parser.parse_args()
    orthoxml_dir = args.orthoxml_dir

    # Add FastOMA utils to the Python path
    add_fastoma_path(args.fastoma_path)

    # Import after updating sys.path
    from OrthoXMLSplitter import OrthoXMLSplitter

    # Change to the directory so that ./splits is created there
    os.chdir(orthoxml_dir)
    print(f"Changed working directory to: {orthoxml_dir}")
    
    # Define the expected filename
    orthoxml_filename = "FastOMA_HOGs.orthoxml"
    # Run the splitter
    splitter = OrthoXMLSplitter(orthoxml_filename, cache_dir="./splits")
    # Extract all HOGs (default behavior: one file per HOG, in ./splits folder)
    splitter()
    logging.info("OrthoXML splitting completed.")

if __name__ == "__main__":
    main()
