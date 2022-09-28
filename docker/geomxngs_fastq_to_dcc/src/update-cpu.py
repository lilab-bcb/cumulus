import argparse
import configparser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--ini", help="INI file", required=True)
    parser.add_argument("--cpu", help="CPU", required=True)
    parser.add_argument("--out", help="Output INI file", required=True)
    args = parser.parse_args()
    ini_path = args.ini
    cpu = args.cpu
    output_path = args.out

    config = configparser.ConfigParser()
    config.optionxform = str  # prevent conversion of keys to lowercase
    config.read(ini_path)
    processing_keys = ["Processing", "Processing_v2"]
    found = False
    for processing_key in processing_keys:
        if processing_key in config:
            config[processing_key]["threads"] = cpu
            found = True
            break
    if not found:
        raise ValueError("Processing section not found")
    with open(output_path, "wt") as out:
        config.write(out)
