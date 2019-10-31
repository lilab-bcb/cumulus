docker build -t pegasus-terra .
docker tag pegasus-terra cumulusprod/pegasus-terra
docker push cumulusprod/pegasus-terra