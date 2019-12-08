docker build -t pegasus-terra-0.15.0 .
docker tag pegasus-terra-0.15.0 cumulusprod/pegasus-terra:0.15.0
docker push cumulusprod/pegasus-terra:0.15.0