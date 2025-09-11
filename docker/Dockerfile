
FROM jupyterlab-datasci-jupyter

RUN apt-get update && apt-get install -y pandoc graphviz
RUN pip3 install pypandoc
RUN pip3 install git+https://github.com/reimerlab/datasci_tools.git

    
EXPOSE 8888