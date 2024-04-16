import random
import math
import matplotlib.pyplot as plt
import numpy as np
class Parameters:
    # Simulation parameters
    lp = 40  # Number of generations in the experiment
    a1 = 1.0  # Initial value of the search space
    b1 = 2.0  # Final value of the search space
    a2 = 3.5  # Initial value of the search space
    b2 = 5.0  # Final value of the search space
    N = 11  # Number of genes in a single chromosome
    pula = 60  # Number of individuals in the population (even number)
    pk = 0.75  # Crossover probability
    pm = 0.01  # Mutation probability
    wykres= [0]*(lp+1)
    counter=0
class BaseNumbers:
    tablica_bazowa = [random.randint(0, 255) for _ in range(100 * Parameters.lp)]
    index_bazowy = 0

    @staticmethod
    def pobierz_bazowa():
        BaseNumbers.index_bazowy += 1
        return BaseNumbers.tablica_bazowa[BaseNumbers.index_bazowy - 1]

    @staticmethod
    def power():
        # Calculate the N-th non-negative power of two
        power = 1
        for _  in range(Parameters.N):
            power *= 2
        return power

class Population:
    def __init__(self):
        self.populacja = [[random.randint(0, 1) for _ in range(Parameters.N)] for _ in range(Parameters.pula)]
        self.tablica_fenotypow = [0] * Parameters.pula
        self.power = BaseNumbers.power()
        self.tablica_dostosowanie = [0] * Parameters.pula

    def losuj_populacje(self):
        for i in range(Parameters.pula):
            for j in range(Parameters.N):
                self.populacja[i][j] = random.randint(0, 1)

    def oblicz_fenotypy(self,a,b):
        for pozycja in range(Parameters.pula):
            self.tablica_fenotypow[pozycja] = a + (b - a) * self.oblicz_fenotyp_chromosomu(pozycja) / self.power

    def oblicz_fenotyp_chromosomu(self, pozycja_chromosomu):
        fenotyp = 0
        rat = 1
        for j in range(Parameters.N):
            fenotyp += self.populacja[pozycja_chromosomu][j] * rat
            rat *= 2
        return fenotyp

    def oblicz_dostosowanie(self,tablica_fenotypowY):
        for i in range(Parameters.pula):
            x = self.tablica_fenotypow[i]
            y= tablica_fenotypowY[i]
            self.tablica_dostosowanie[i] = 1-math.log(x**2+math.cos(y),10)

    def dostosowanie_normalizacja(self):
        min_val = min(self.tablica_dostosowanie)
        max_val = max(self.tablica_dostosowanie)
        offset = (max_val - min_val) / (Parameters.N - 1) - min_val
        for i in range(Parameters.pula):
            self.tablica_dostosowanie[i] += offset

    def ruletka(self):
        self.dostosowanie_normalizacja()
        suma_dostosowanie = sum(self.tablica_dostosowanie)

        tablica_NI = [dostosowanie / suma_dostosowanie * self.power for dostosowanie in self.tablica_dostosowanie]
        losowe = [random.randint(0, self.power - 1) for _ in range(Parameters.pula)]
        ruletka = [0] * Parameters.pula

        pozycja = 0
        for i in range(Parameters.pula):
            pozycja += tablica_NI[i]
            ruletka[i] = pozycja

        nowe_pokolenie = [[0] * Parameters.N for _ in range(Parameters.pula)]

        for i in range(Parameters.pula):
            j = 0
            while losowe[i] > ruletka[j]:
                j += 1
            nowe_pokolenie[i] = self.populacja[j][:]

        self.populacja = nowe_pokolenie

    def krzyzowanie(self):
        liczba_par = Parameters.pula // 2
        losowe_pary = [random.randint(0, 99) for _ in range(liczba_par)]
        losowe_miejsca = [random.randint(0, Parameters.N - 2) for _ in range(liczba_par)]

        pierwszy = 0
        for para in range(liczba_par):
            if losowe_pary[para] < Parameters.pk * 100:
                for i in range(losowe_miejsca[para], Parameters.N):
                    self.populacja[pierwszy][i], self.populacja[pierwszy + 1][i] = self.populacja[pierwszy + 1][i], self.populacja[pierwszy][i]
            pierwszy += 2

    def mutacje(self):
        for i in range(Parameters.pula):
            if random.random() < Parameters.pm:
                miejsce_mutacji = random.randint(0, Parameters.N - 1)
                self.populacja[i][miejsce_mutacji] = 1 - self.populacja[i][miejsce_mutacji]

    def pokaz_dostosowanie_srednie(self):
        srednia = sum(self.tablica_dostosowanie) / Parameters.pula
        Parameters.wykres[Parameters.counter]=srednia
        Parameters.counter+=1
        print(f"{srednia:.2f}")

def main():
    nr_pokolenia = 0
    liczby_bazowe = BaseNumbers()
    populacja1 = Population()
    populacja2 = Population()

    populacja1.losuj_populacje()
    populacja1.oblicz_fenotypy(Parameters.a1,Parameters.b1)

    populacja2.losuj_populacje()
    populacja2.oblicz_fenotypy(Parameters.a2,Parameters.b2)

    populacja1.oblicz_dostosowanie(populacja2.tablica_fenotypow)
    populacja2.tablica_dostosowanie=populacja1.tablica_dostosowanie


    print("Nr pokolenia   Średnia wartość funkcji dostosowania")
    print(f"{nr_pokolenia:3}", end="          ")
    populacja1.pokaz_dostosowanie_srednie()
    while nr_pokolenia < Parameters.lp:
        nr_pokolenia += 1
        populacja1.ruletka()
        populacja1.krzyzowanie()
        populacja1.mutacje()
        populacja1.oblicz_fenotypy(Parameters.a1,Parameters.b1)




        populacja2.ruletka()
        populacja2.krzyzowanie()
        populacja2.mutacje()
        populacja2.oblicz_fenotypy(Parameters.a2,Parameters.b2)

        populacja1.oblicz_dostosowanie(populacja2.tablica_fenotypow)
        populacja2.tablica_dostosowanie = populacja1.tablica_dostosowanie
        print(f"{nr_pokolenia:3}", end="          ")
        populacja1.pokaz_dostosowanie_srednie()


    numbers = np.arange(1, 42)
    plt.plot(numbers, Parameters.wykres,c = "b")

    # Set the title
    plt.title("Algorytm genetyczny dla funkcji   1-Log[10,x^2+cos(y)]  ")
    # Set the y-axis label
    plt.ylabel('wartosc funkcji')
    # Set the x-axis label
    plt.xlabel('Numer pokolenia')
    plt.show()
if __name__ == "__main__":
    main()



