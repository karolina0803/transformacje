### PROJEKT 1: TRANSFORMACJE

Program służy do transformacji współrzędnych.
## Spis treści
* [Dostępne transformacje](#dostępne-transformacje)
* [Dostępne elipsoidy](#dostępne-elipsoidy)
* [Wymagania](#wymagania)
* [Opis programu](#opis-programu)
* [Wynik końcowy](#wynik-końcowy)
* [Przykład użycia](#przykładowe-transformacje)
* [Znane błędy](#znane-błędy)

#### Dostępne transformacje:
```
XYZ -> PLH (Zamiana współrzędnych geocentrycznych na elipsoidalne)
PLH -> XYZ (Zamiana współrzędnych elipsoidalnych na geocentryczne)
XYZ -> NEU (Zamiana współrzędnych geocentrycznych na współrzędne w układzie topocentrycznym)
PLH -> PL2000 (Zamiana współrzędnych elipsoidalnych na współrzędne w układzie PL-2000)
PLH -> PL1992 (Zamiana współrzędnych elipsoidalnych na współrzędne w układzie PL-1992)
XYZkras -> XYZgrs80 (Zamiana współrzędnych z elipsoidy Krasowskiego na GRS80)
XYZgrs80 -> PL2000 (Zamiana współrzędnych na elipsoidzie GRS80 na współrzędne w układzie PL2000) 
XYZgrs80 -> PL1992 (Zamiana współrzędnych na elipsoidzie GRS80 na współrzędne w układzie PL1992)
```
#### Dostępne elipsoidy

cokolwiek
```
GRS80
WGS84
Elipsoida Krasowskiego
```
#### Wymagania
```
Python 3.11
Biblioteka numpy
Biblioteka math
System operacyjny Windows 10 lub 11
```
### Opis programu
Wybór transformacji odbywa się po wpisaniu jednej z następujących flag:
```
--xyz2plh (transformacja XYZ->PLH)
--plh2xyz (transformacja PLH->XYZ)
--xyz2neu (transformacja XYZ->NEU)
--plh2pl2000 (transformacja PLH->PL-2000)
--plh2pl1992 (transformacja PLH->PL1992)
--XYZkras2XYZgrs80 
--XYZgrs802pl2000  
--XYZgrs802pl1992
```
Wybór elipsoidy możliwy jest poprzez wpisanie tylko jednej z następujących nazw:
```
grs80
wgs84
elipsoida_krasowskiego
```
Do programu należy przekazać plik tekstowy ze współrzędnymi. Przykładowy format dla współrzędnych XYZ: 
XXX.XXX,YYY.YYY,ZZZ.ZZZ

Plik ze współrzędnymi musi znajdować się w tym samym folderze co program do transformacji.

Jeśli plik ze współrzędnymi zawiera sekcję nagłówka, to w programie należy wpisać jako argument liczbę linijek tego nagłówka.

Argumenty należy podawać w kolejności:
```
1.python
2.nazwa programu
3.liczba linijek nagłówka
4.model elipsoidy
5.rodzaj transformacji
6.nazwa pliku ze współrzędnymi

```
##### UWAGA
W celu wykonania transformacji XYZ->NEU należy ręcznie wpisać współrzędne geocentryczne stanowiska po argumencie z rodzajem transformacji w kolejności XXX.XXX YYY.YYY ZZZ.ZZZ, a dopiero potem wpisać nazwę pliku ze współrzędnymi.

#### Wynik końcowy
Po poprawnym wpisaniu wszystkich argumentów wyświetli się komunikat:
```
 korzystasz z modelu [model elipsoidy]

transformacja przebigla poprawnie
 plik wynikowy dostepny pod nazwa 'result_[transformacja].txt'
```
[model elipsoidy] to wybrany model do transformacji

[transformacja] to nazwa transformacji, która została wykonana

#### Przykład użycia
```
python transformacje_1.py 4 grs80 --xyz2plh wsp.txt
```
lub w przypadku XYZ->NEU:
```
python transformacje_1.py 4 grs80 --xyz2neu 100.000 100.000 100.000 wsp.txt
```
#### Przykładowe transformacje

XYZ -> BLH dla poniższych danych:
```
3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
```
Wynik wyświetli się w następujący sposób w pliku o nazwie result_XYZ_philamh.txt:

```
    phi [deg]      lam [deg]         h [m]
-----------------------------------------------
   52.097272   ,   21.031533   ,    141.399    
   52.097272   ,   21.031533   ,    141.400 
  
```
XYZ -> NEU dla danych:
```
3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
```
Współrzędne stanowiska to: 100.000 100.000 100.000

Wynik wyśwetli się w następujący sposób w pliku o nazwie result_XYZ_NEU.txt:
```
     N [m]          E [m]           U [m]
--------------------------------------------
 5017888.028  ,  -1595082.221 ,  3576003.602  
 5017888.025  ,  -1595082.235 ,  3576003.602  
```
PLH->PL2000 dla danych:
```
   52.097272   ,   21.031533   ,    141.399    
   52.097272   ,   21.031533   ,    141.400 
```
Wynik wyśwetli się w następujący sposób w pliku o nazwie result_philamh_pl2000.txt:
```
     X [m]             Y [m]    
---------------------------------
  5773722.697   ,   7502160.760  
  5773722.697   ,   7502160.760 
```

#### Znane błędy

1. Przy transformacji XYZ->PLH, a następnie PLH->XYZ otrzymujemy wyniki XYZ różniące się o 2 cm względem pliku wejściowego



