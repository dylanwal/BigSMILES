# Testing

## Pytest

Test function of code against known answers. 

Run `main.py`.

Run with coverage report `pytest --cov --cov-report=html:coverage_re`


## performance

Test CPU time and memory usage with version updates.

Run `performance_analysis\main.py`.

Results are appended on to 'performance.csv'.



# Building

`python -m build`

`twine upload dist/*`
