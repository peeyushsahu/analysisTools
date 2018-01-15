FROM python:3.4
RUN mkdir /python_pkg
RUN mkdir /python_pkg/results
COPY /python_pkg/ /python_pkg/
CMD ["python", "python_pkg/test_script.py"]

