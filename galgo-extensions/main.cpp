#include "Galgo.hpp"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template <typename T>
vector<T> objective(const vector<T>& x) {
    // The count of satisfied conditions
    int conditions = 0;
    int i;

    // 1. Ron jí zmrzlinu
    if (static_cast<int>(x[0]) == 1) conditions++;

    // 2. Joe nemá rád žvýkačky
    if (static_cast<int>(x[1]) != 0) conditions++;

    // 3. Samovi je 14
    if (static_cast<int>(x[12]) == 14) conditions++;

    // 4. Sam není na horské dráze
    if (static_cast<int>(x[7]) != 0) conditions++;

    // 5. Chlapci na zvonkové dráze je 15
    // hledám, kdo je na zvonkové dráze
    for (i=5;i<10;i++)
      if (static_cast<int>(x[i])==1) break;
    if (i!=10 && static_cast<int>(x[i+5])==15 ) conditions++;

    // 6. Len není ve strašidelném zámku
    if (static_cast<int>(x[8]) != 2) conditions++;

    // 7. Don je na kolotoči
    if (static_cast<int>(x[9]) == 3) conditions++;

    // 8. Chlapci, který nemá rád zmrzlinu, je 13 let
    // hledam, komu je 13
    for (i=10;i<15;i++)
      if (static_cast<int>(x[i])==13) break;
    if (i!=15 && static_cast<int>(x[i-10])!=1) conditions++;

    // 9. Chlapec ve strašidelném zámku jí párek v rohlíku
    // hledam, kdo je ve strašidelném zámku
    for (i=5;i<10;i++)
      if (static_cast<int>(x[i])==2) break;
    if (i!=10 && static_cast<int>(x[i-5])==2) conditions++;

    // 10. Joe je jedenáctiletý
    if (static_cast<int>(x[11]) == 11) conditions++;

    // 11. Joe jí hranolky
    if (static_cast<int>(x[1]) == 3) conditions++;

    // 12. Joe je na lochnesce
    if (static_cast<int>(x[6]) == 4) conditions++;

    // 13. Donovi je 12 let
    if (static_cast<int>(x[14]) == 12) conditions++;

    // 14. Don si pochutnává na cukrové vatě
    if (static_cast<int>(x[4]) == 4) conditions++;

    T objectiveValue = static_cast<T>(conditions);

    // Vytvoření vektoru s hodnotou objektivní funkce
    vector<T> result = { objectiveValue };

    return result;
}

vector<double> MyConstraint(const vector<double>& x)
{
   vector<double> cst;
   int index[3]={0,5,10};
   int val;
   int i;

   for(i=0;i<3;i++)
   {
     val=-1;
     if (static_cast<int>(x[index[i]+0])==static_cast<int>(x[index[i]+1])) val=14;
     if (static_cast<int>(x[index[i]+0])==static_cast<int>(x[index[i]+2])) val=14;
     if (static_cast<int>(x[index[i]+0])==static_cast<int>(x[index[i]+3])) val=14;
     if (static_cast<int>(x[index[i]+0])==static_cast<int>(x[index[i]+4])) val=14;
     if (static_cast<int>(x[index[i]+1])==static_cast<int>(x[index[i]+2])) val=14;
     if (static_cast<int>(x[index[i]+1])==static_cast<int>(x[index[i]+3])) val=14;
     if (static_cast<int>(x[index[i]+1])==static_cast<int>(x[index[i]+4])) val=14;
     if (static_cast<int>(x[index[i]+2])==static_cast<int>(x[index[i]+3])) val=14;
     if (static_cast<int>(x[index[i]+2])==static_cast<int>(x[index[i]+4])) val=14;
     if (static_cast<int>(x[index[i]+3])==static_cast<int>(x[index[i]+4])) val=14;
     cst.push_back(static_cast<double>(val));
   }
  return cst;
}

template <typename T>
void ProblemCrossover(const galgo::Population<T>& x, galgo::CHR<T>& chr1, galgo::CHR<T>& chr2)
{
   // choosing randomly 2 chromosomes from mating population
   int idx1 = galgo::uniform<int>(0, x.matsize());
   int idx2 = galgo::uniform<int>(0, x.matsize());
   // choosing randomly a position for cross-over between two boundaries
   int pos = galgo::uniform<int>(1, 3);
   pos = pos*5*3-1;
   // transmitting portion of bits to new chromosomes
   chr1->setPortion(*x[idx1], 0, pos);
   chr2->setPortion(*x[idx2], 0, pos);
   chr1->setPortion(*x[idx2], pos + 1);
   chr2->setPortion(*x[idx1], pos + 1);
}


template <typename T>
void ProblemMutation(galgo::CHR<T>& chr)
{
   T mutrate = chr->mutrate();

   if (mutrate == 0.0) return;

   // swap mutation only within group of 5 elements
   for (int i = 0; i < chr->nbgene(); ++i) {
      // Apply mutation with probability mutrate
      if (galgo::proba(galgo::rng) <= mutrate) {
         // Generate random indices for swapping
         int low = (i / 5)*5; // lower bound of swap boundaries
         int high = low + 5; // higher bound of swap boundaries
         int index1 = galgo::uniform<int>(low, high);
         int index2 = galgo::uniform<int>(low, high);
         //Swap the genes at index1 and index2
         double gen1=chr->getGene(index1);
         double gen2=chr->getGene(index2);
         chr->setGeneValue(index1,gen2);
         chr->setGeneValue(index2,gen1);
      }
   }
}

//template <typename T>
//void MyMutation2(galgo::CHR<T>& chr)
//{
//   T mutrate = chr->mutrate();
//
//   if (mutrate == 0.0) return;
//
//   // swap mutation only within group of 5 elements
//   if (galgo::proba(galgo::rng) <= mutrate) {
//         // Generate random indices for swapping
//         int index1 = galgo::uniform<int>(0, chr->nbgene());
//         int group = index1/5;
//         int index2 = galgo::uniform<int>(0,5);
//         index2=5*group+index2;
//         // Swap the genes at index1 and index2
//         std::swap(chr->param[index1], chr->param[index2]);
//      }
//}
const char *convert_food(double ind)
{
    int pom = static_cast<int>(ind);
  switch (pom) {
    case 0: return "Zvykacka";
    case 1: return "Zmrzlina";
    case 2: return "Parek v rohliku";
    case 3: return "Hranolky";
    case 4: return "Cukrova vata";
  }
  return "";
}

const char *convert_attraction(double ind)
{
    int pom = static_cast<int>(ind);
  switch (pom) {
    case 0: return "Horska draha";
    case 1: return "Zvonkova draha";
    case 2: return "Strasidelny zamek";
    case 3: return "Kolotoc";
    case 4: return "Lochneska";
  }
  return "";
}

int convert_age(double ind)
{
    int pom = static_cast<int>(ind);
  return pom;
}

const char *convert_name(int ind)
{
  switch (ind) {
    case 0: return "Ron";
    case 1: return "Joe";
    case 2: return "Sam";
    case 3: return "Len";
    case 4: return "Don";
 }
 return "";
}

int main() {
    // Initialize parameters lower and upper bounds
    // Define parameter bounds for food
    galgo::Parameter<double,3> RonFood({0, 4.9,1.2}); // 0: Žvýkačka, 1: Zmrzlina, 2: Párek v rohlíku, 3: Hranolky, 4: Cukrová vata
    galgo::Parameter<double,3> JoeFood({0, 4.9,3.5}); // 0: Žvýkačka, 1: Zmrzlina, 2: Párek v rohlíku, 3: Hranolky, 4: Cukrová vata
    galgo::Parameter<double,3> SamFood({0, 4.9,0.5}); // 0: Žvýkačka, 1: Zmrzlina, 2: Párek v rohlíku, 3: Hranolky, 4: Cukrová vata
    galgo::Parameter<double,3> LenFood({0, 4.9,2.8}); // 0: Žvýkačka, 1: Zmrzlina, 2: Párek v rohlíku, 3: Hranolky, 4: Cukrová vata
    galgo::Parameter<double,3> DonFood({0, 4.9,4.2}); // 0: Žvýkačka, 1: Zmrzlina, 2: Párek v rohlíku, 3: Hranolky, 4: Cukrová vata
    // Define parameter bounds for attractions
    galgo::Parameter<double,3> RonAttraction({0, 4.9,0.3}); // 0: Horská dráha, 1: Zvonková dráha, 2: Strašidelný zámek, 3: Kolotoč, 4: Lochneska
    galgo::Parameter<double,3> JoeAttraction({0, 4.9,4.3}); // 0: Horská dráha, 1: Zvonková dráha, 2: Strašidelný zámek, 3: Kolotoč, 4: Lochneska
    galgo::Parameter<double,3> SamAttraction({0, 4.9,2.7}); // 0: Horská dráha, 1: Zvonková dráha, 2: Strašidelný zámek, 3: Kolotoč, 4: Lochneska
    galgo::Parameter<double,3> LenAttraction({0, 4.9,1.6}); // 0: Horská dráha, 1: Zvonková dráha, 2: Strašidelný zámek, 3: Kolotoč, 4: Lochneska
    galgo::Parameter<double,3> DonAttraction({0, 4.9,3.2}); // 0: Horská dráha, 1: Zvonková dráha, 2: Strašidelný zámek, 3: Kolotoč, 4: Lochneska
    // Define parameter bounds for age
    galgo::Parameter<double,3> RonAge({11, 15.9,15.6}); // Věk od 11 do 15 let
    galgo::Parameter<double,3> JoeAge({11, 15.9,11.6}); // Věk od 11 do 15 let
    galgo::Parameter<double,3> SamAge({11, 15.9,14.7}); // Věk od 11 do 15 let
    galgo::Parameter<double,3> LenAge({11, 15.9,13.8}); // Věk od 11 do 15 let
    galgo::Parameter<double,3> DonAge({11, 15.9,12.8}); // Věk od 11 do 15 let

    // Initialize genetic algorithm with defined parameters
    galgo::GeneticAlgorithm<double> ga(objective,800,500,true,0,RonFood, JoeFood, SamFood, LenFood,DonFood,RonAttraction,JoeAttraction,SamAttraction,LenAttraction,DonAttraction,RonAge,JoeAge,SamAge,LenAge,DonAge);

    ga.CrossOver = ProblemCrossover;
    ga.Mutation=ProblemMutation;
    // The following comment were used for testing with different (swap) mutation
    //ga.Mutation=SwapMutation;
    ga.covrate = 0.8;
    ga.mutrate = 0.1;
    ga.elitpop=100;
    ga.Selection=TNT;
    ga.tntsize=2;
    ga.Constraint=MyConstraint;
    ga.FitnessBestTarget = 14.0;
    // Run genetic algorithm
    ga.run();

	std::cout << "Nalezene reseni:\n";

	int i;
    for(i = 0; i < 5; i++)   std::cout << convert_name(i) << ": " << convert_food(ga.pop(0)->getParam()[i]) << std::endl;
    for(i = 5; i < 10; i++)  std::cout << convert_name(i - 5) << ": " << convert_attraction(ga.pop(0)->getParam()[i]) << std::endl;;
    for(i = 10; i < 15; i++) std::cout << convert_name(i - 10) << ": " << convert_age(ga.pop(0)->getParam()[i]) << std::endl;;;

    return 0;
}
