FROM python:3.4

RUN mkdir /GCAM
RUN mkdir /resource
RUN mkdir /results

COPY analysisTools/ /GCAM

WORKDIR /GCAM

RUN pip install -r requirement.txt

WORKDIR /GCAM/ 

ENTRYPOINT ["/bin/sh", "run.sh"]
CMD []
