# run the follwoing command and use the localhost:8888 on your browser
#
ssh -NL 1234:localhost:8888 centos@10.10.1.22

# Run the following and then http://localhost:9000
ssh -i ~/ec2-keys/cent-os-keys.cer  -NL 9000:localhost:8888 ubuntu@10.10.1.29
