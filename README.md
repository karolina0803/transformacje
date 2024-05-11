### PROJEKT 1: TRANSFORMACJE

Program służy do transformacji współrzędnych.
#### Dostępne transformacje:
```
XYZ -> PLH (Zamiana współrzędnych geocentrycznych na elipsoidalne)
PLH -> XYZ (Zamiana współrzędnych elipsoidalnych na geocentryczne)
XYZ -> NEU (Zamiana współrzędnych geocentrycznych na współrzędne w układzie topocentrycznym)
PLH -> PL2000 (Zamiana współrzędnych elipsoidalnych na współrzędne w układzie PL-2000)
PLH -> PL1992 (Zamiana współrzędnych elipsoidalnych na współrzędne w układzie PL-1992)
???? XYZkras -> XYZgrs80 (Zamiana współrzędnych z elipsoidy Krasowskiego na GRS80)
???? XYZgrs80 -> PL2000 (Zamiana współrzędnych na elipsoidzie GRS80 na współrzędne w układzie PL2000) 
???? XYZgrs80 -> PL1992 (Zamiana współrzędnych na elipsoidzie GRS80 na współrzędne w układzie PL1992)
```
#### Dostępne elipsoidy
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
--XYZkras2XYZgrs80 ????
--XYZgrs802pl2000  ????
--XYZgrs802pl1992' ????
```
Wybór elipsoidy możliwy jest poprzez wpisanie tylko jednej z następujących nazw:
```
grs80
wgs84
elipsoida_krasowskiego
```
Do programu należy przekazać plik tekstowy ze współrzędnymi. Przykładowy format dla współrzędnych XYZ: 
XXX.XXX,YYY.YYY,ZZZ.ZZZ

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

````
python transformacje_1.py 4 grs80 --xyz2plh wsp.txt
```
lub w przypadku XYZ->NEU:

