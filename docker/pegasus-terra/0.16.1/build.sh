docker build -t pegasus-terra-0.16.1 .
docker tag pegasus-terra-0.16.1 cumulusprod/pegasus-terra:0.16.1
docker push cumulusprod/pegasus-terra:0.16.1