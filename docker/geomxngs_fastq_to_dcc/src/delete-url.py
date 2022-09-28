import argparse
from subprocess import check_call


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--backend", required=True)
    parser.add_argument("--url", required=True)

    args = parser.parse_args()
    backend = args.backend
    url = args.url
    if not url.endswith("/"):
        url += "/"
    try:
        call_args = ["strato", "rm", "--backend", backend, "-m", "-r", url]
        check_call(call_args)
        print("Deleted {}".format(url))
    except:
        print("Failed to delete {}.".format(url))
