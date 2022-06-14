import sys

if __name__ == "__main__":
    backend = sys.argv[1]
    uri_name = sys.argv[2]

    if backend == 'gcp':
        assert uri_name.startswith('gs://'), f"Your output directory '{uri_name}' does not match {backend} backend"
    elif backend == 'aws':
        assert uri_name.startswith('s3://'), f"Your output directory '{uri_name}' does not match {backend} backend"
    else:
        assert (not uri_name.startswith('s3://')) and (not uri_name.startswith('gs://')), f"Your output directory '{uri_name}' does not match {backend} backend"
