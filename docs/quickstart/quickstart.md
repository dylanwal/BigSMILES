This part of the project documentation focuses on a **learning-oriented** approach. 
- Help newcomers with getting started


suggested shorting of `import bigsmiles` --> `import bigsmiles as bs`


```python
import bigsmiles

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bigsmiles.BigSMILES(polymer_string)
print(bigsmiles)  # "CC{[>][<]CC(C)[>][<]}CC(C)=C"
```