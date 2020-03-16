#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include <vector>
#include <cmath>


using namespace std;

#define PI 3.14159265

float Dixon_Price(vector<float> x, int size)
{
    float result = pow(x[0] - 1, 2);

    for(int i = 1; i < size; ++i)
    {
        result += i * pow (( 2 * pow(x[i], 2) - x[i-1]), 2);
    }

    return result;
}

float Rosenbrock(vector<float> x, int size)
{
    float result = 0;

    for(int i = 0; i < size - 1; ++i)
    {
        result += 100 * pow (x[i+1] - pow(x[i], 2), 2) + pow(x[i]-1, 2);
    }

    return result;
}

float sphere(vector<float> x, int size)
{
    float result = 0;

    for(int i = 0; i < size; i++)
    {
        result += pow(x[i],2);
    }
    return result;
}

float Rastrigin(vector<float> x, int size)
{
    float result = 10 * size;

    for(int i = 0; i < size; ++i)
    {
        result += pow(x[i],2) - 10 * cos(2*PI*x[i]);
    }
    return result;
}

void afisare(vector<int> vect, int D, int L){
    // afisare biti solutie
    for (int i = 0; i < D*L; i++) {
        if(i % L == 0 && i > 0) cout<< "\n" ;
        cout<< vect[i]<< " ";
    }
    cout<<endl;
}

float random01()
{
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

float convert_binary(int x_binar[], int L, float a, float b)
{
    float x = 0;
    /*
    for(int i = 0; i < L; ++i)
    {
        cout<<x_binar[i];
    }
    cout<<endl;
     */
    for(int i = 0; i < L; ++i)
    {
        x = 2 * x;
        x += x_binar[i];
    }
    //cout<<endl;
    x = x / pow(2,L);
    x *= (b-a);
    x += a;
    return x;
}

vector<float> deco(vector<int> x_binar, int L, int D, float a, float b)
{
    vector<float> x;
    int vec[L];
    for(int i = 0, j = 0; i < L*D; ++i, ++j){
        if(i % L == 0 && i > 0){
            j = 0;
            x.push_back(convert_binary(vec, L, a, b));
        }
        vec[j] = x_binar.at(i);
    }
    x.push_back(convert_binary(vec, L, a, b));
    return x;
}

float eval(vector<int> binar, int D, int L, int funct)
{
    vector<float> x_arg;
    float result = 0;
    if(funct == 0)
    {
        float a = -10, b = 10;
        x_arg = deco(binar, L, D, a, b);
        result = Dixon_Price(x_arg, D);
    }
    if(funct == 1)
    {
        float a = -5, b = 10;
        x_arg = deco(binar, L, D, a, b);
        result = Rosenbrock(x_arg, D);
    }
    if(funct == 2)
    {
        float a = -5.12, b = 5.12;
        x_arg = deco(binar, L, D, a, b);
        result = sphere(x_arg, D);
    }
    if(funct == 3)
    {
        float a = -5.12, b = 5.12;
        x_arg = deco(binar, L, D, a, b);
        result = Rastrigin(x_arg, D);
    }
    return result;
}

vector<vector<int> > generate_population(int D, int L, int pop_size)
{
    vector<vector<int> > population;
    
    for(int p = 0; p < pop_size; ++p)
    {
        //randomly generate one individual
        vector<int> individual;
        for (int i = 0; i < D * L; ++i) {
                individual.push_back(rand() % 2);
        }
        
        //add the individual to the population
        population.push_back(individual);
    }
    return population;
}

vector<float> fitness(vector<vector<int> > population, int D, int L, int funct)
{
    float maxF = -1;
    vector<int> individual;
    vector<float> x_arg;
    vector<float> fitnessP;
    vector<float> functVal;
    
    //calculate maximum function value of population
    for(int i = 0; i < population.size(); ++i)
    {
        individual = population.at(i);
        if(funct == 0)
        {
            float a = -10;
            float b = 10;
            
            //decode the individual
            x_arg = deco(individual, L, D, a, b);
            //calculate value of function on that individual
            float eval_funct = Dixon_Price(x_arg, D);
            //memorate value of function
            functVal.push_back(eval_funct);
            //compare maximum with value found
            if(maxF < eval_funct)
            {
                //change maximum
                maxF = eval_funct;
            }
        }
        if(funct == 1)
        {
            float a = -5;
            float b = 10;
            
            //decode the individual
            x_arg = deco(individual, L, D, a, b);
            //calculate value of function on that individual
            float eval_funct = Rosenbrock(x_arg, D);
            //memorate value of function
            functVal.push_back(eval_funct);
            //compare maximum with value found
            if(maxF < eval_funct)
            {
                //change maximum
                maxF = eval_funct;
            }
        }
        if(funct == 2)
        {
            float a = -5.12;
            float b = 5.12;
            
            //decode the individual
            x_arg = deco(individual, L, D, a, b);
            //calculate value of function on that individual
            float eval_funct = sphere(x_arg, D);
            //memorate value of function
            functVal.push_back(eval_funct);
            //compare maximum with value found
            if(maxF < eval_funct)
            {
                //change maximum
                maxF = eval_funct;
            }
        }
        if(funct == 3)
        {
            float a = -5.12;
            float b = 5.12;
            
            //decode the individual
            x_arg = deco(individual, L, D, a, b);
            //calculate value of function on that individual
            float eval_funct = Rastrigin(x_arg, D);
            //memorate value of function
            functVal.push_back(eval_funct);
            //compare maximum with value found
            if(maxF < eval_funct)
            {
                //change maximum
                maxF = eval_funct;
            }
        }
        x_arg.clear();
        individual.clear();
    }
    
    //float epsilon = pow(10,-9);
    //calculate fitness for each individual
    for(int i = 0; i < population.size(); ++i)
    {
        //fitnessP.push_back(1/functVal.at(i) + epsilon);
        fitnessP.push_back(1.1 * maxF - functVal.at(i));
    }

    return fitnessP;
}

vector<float> calculate_sf(vector<float> fitnessP)
{
    vector<float> sf;
    
    sf.push_back(fitnessP.at(0));
    for(int i = 1; i < fitnessP.size(); ++i)
    {
        //compute probability of selecting individual i
        sf.push_back(sf.at(i-1) + fitnessP.at(i));
    }
    return sf;
}

int selectF(vector<float> sf)
{
    float pos =  random01() * sf.at(sf.size()-1);
    
    int i;
    for(i = 0; i < sf.size(); ++i)
    {
        if(pos <= sf.at(i))
            return i;
    }
    return i;
}

vector<vector<int> > selection(vector<vector<int> > population, vector<float> fitnessP)
{
    vector<vector<int> > new_pop;
    vector<float> sf = calculate_sf(fitnessP);
    vector<int> index;
    for(int i = 0; i < population.size(); ++i)
    {
        int j = selectF(sf);
        new_pop.push_back(population.at(j));
    }
    return new_pop;
}

void mutation(vector<vector<int> >& population, int generation)
{
    vector<int> individual;
    float pm;
    
    if(generation < 600){
        pm = 0.01;
    }
    else{
        pm = 0.02;
    }

    
    for(int i = 0; i < population.size(); ++i)
    {
        individual = population.at(i);
        for(int j = 0; j < population.at(i).size(); ++j)
        {
            if(random01() < pm)
            {
                individual.at(j) = !individual.at(j);
            }
        }
        
        population.at(i) = individual;
        individual.clear();
    }
}
void  crossover(vector<vector<int> >& population, int L, int D, int generation)
{
    float pc;
    if(generation < 600){
        pc = 0.3;
    }else{
        pc = 0.25;
    }
    vector<int> pos;
    vector<int> individual1;
    vector<int> individual2;
    vector<int> copyind1;
    vector<int> copyind2;
    
    //choose 2 random individuals
    for(int i = 0; i < population.size(); ++i)
    {
        if(random01() < pc)
        {
            pos.push_back(i);
        }
    }
    //work on prob
    if(pos.size() % 2 == 1)
    {
        pos.pop_back();
    }
    //do cross over
    for(int i = 0; i < pos.size(); i +=2)
    {
        individual1 = copyind1 = population.at(pos.at(i));
        individual2 = copyind2 = population.at(pos.at(i+1));
        
        int slice = rand() % (D * L -2) + 1;
        
        for(int j = slice; j < D * L; ++j)
        {
            copyind1.at(j) = individual2.at(j);
            copyind2.at(j) = individual1.at(j);
        }
        //get new values
        population.at(pos.at(i)) = copyind1;
        population.at(pos.at(i+1)) = copyind2;
        
        individual1.clear();
        individual2.clear();
        copyind1.clear();
        copyind2.clear();
    }
}

int get_position_best(vector<float> fit)
{
    float maxValue = -1;
    int pos = 0;
    
    for(int i = 0; i < fit.size(); ++i)
    {
        if(fit.at(i) > maxValue)
        {
            maxValue = fit.at(i);
            pos = i;
        }
    }
    return pos;
}

float get_best_value(vector<vector<int> > population, vector<float> fit, int funct, int D, int L)
{
    vector<float> x_arg;
    vector<int> individual;
    int posmin;
    posmin = get_position_best(fit);
    individual = population.at(posmin);
    if(funct == 0)
    {
        float a = -10, b = 10;
        //decode the individual
        x_arg = deco(individual, L, D, a, b);
        //calculate value of function on that individual
        float eval_funct = Dixon_Price(x_arg, D);
        return eval_funct;
    }
    if(funct == 1)
    {
        float a = -5, b = 10;
        //decode the individual
        x_arg = deco(individual, L, D, a, b);
        //calculate value of function on that individual
        float eval_funct = Rosenbrock(x_arg, D);
        return eval_funct;
    }
    if(funct == 2)
    {
        float a = -5.12, b = 5.12;
        //decode the individual
        x_arg = deco(individual, L, D, a, b);
        //calculate value of function on that individual
        float eval_funct = sphere(x_arg, D);
        return eval_funct;
    }
    if(funct == 3)
    {
        float a = -5.12, b = 5.12;
        //decode the individual
        x_arg = deco(individual, L, D, a, b);
        //calculate value of function on that individual
        float eval_funct = Rastrigin(x_arg, D);
        return eval_funct;
    }
    return 0;
}

vector<int> neighbor(vector<int> x_binar, int pos)
{
    vector<int> candidat;
    candidat = x_binar;
    candidat[pos] = !candidat[pos];
    return candidat;
}



vector<int> Improve(vector<int> x_binar, int L, int D, int funct)
{
    float min_gasit;
    vector<float> x_arg;
    float fc;
    vector<int> candidat;
    vector<float> arg_cand;
    vector<int> best;
    int optiune = 1;

    min_gasit = eval(x_binar, D, L, funct);
    
    for(int i = 0; i < L*D; ++i)
    {
        candidat = neighbor(x_binar, i);
        
        //calculam functia in candidat
        fc = eval(candidat, D, L, funct);
        
        if(fc < min_gasit && optiune == 0)
        {
            min_gasit = fc;
            best = candidat;
            break;
        }
        if(fc < min_gasit && optiune == 1)
        {
            min_gasit = fc;
            best = candidat;
        }
    }
    if(best.empty())
    {
        return x_binar;
    }
    //cout<<min_gasit;
    return best;
}



float Hill_Climbing(int L, int D, vector<int> best, int funct)
{
    
    float candidat_value;
    float best_value;
    float new_sol_value;
    vector<int> candidat;
    vector<int> new_sol;
    
    //evaluam valoarea functiei
    best_value = eval(best, D, L, funct);
    cout<<best_value<<endl;
    
    int t = 100;
    while(t)
    {
        int local = 1;
        
        if(t == 100){
            candidat = best;
            candidat_value = best_value;
        }else{
            //select candidate
        for (int i = 0; i < D * L; i++) {
            candidat.push_back(rand() % 2);
        }
        // eval candidat
        candidat_value = eval(candidat, D, L, funct);
        
        }

        while(local){
            
            new_sol = Improve(candidat, L, D, funct);
            new_sol_value = eval(new_sol, D, L, funct);
    
            if(new_sol_value < candidat_value)
            {
                candidat.clear();
                candidat = new_sol;
                candidat_value = new_sol_value;
            }
            else
            {
                local = 0;
            }
        }
        
        if(candidat_value < best_value)
        {
            best.clear();
            best = candidat;
            best_value = candidat_value;
        }
        t = 0;
        candidat.clear();
    }
    
    return best_value;
}

float genetic_algorithm(int D, int L, int generations, int funct, int pop_size, ofstream &file)
{
    
    float minFunc;
    vector<int> individual;
    int hill = 1;
    //generate intial population
    vector<vector<int> > population = generate_population(D, L, pop_size);
    //calculate fitness
    vector<float> fit = fitness(population, D, L, funct);
    //get best out of this population
    minFunc = get_best_value(population, fit, funct, D, L);
    
    if(hill == 1){
        individual.clear();
        int posmin;
        posmin = get_position_best(fit);
        individual = population.at(posmin);
    }
    
    while(generations)
    {
        //alter
        mutation(population, generations);
        crossover(population, L, D, generations);
        
        //calcultate new fitness
        fit.clear();
        fit = fitness(population, D, L, funct);
        
        //selection
        vector<vector<int> > new_population = selection(population, fit);
        population.clear();
        population = new_population;
        
        //calcultate new fitness
        fit.clear();
        fit = fitness(population, D, L, funct);
        
        //evaluate
        float minGeneration = get_best_value(population, fit, funct, D, L);
        //cout <<generations << " " << minGeneration << endl;
        
        if(minFunc > minGeneration)
        {
            minFunc = minGeneration;
            if(hill == 1){
                individual.clear();
                int posmin;
                posmin = get_position_best(fit);
                individual = population.at(posmin);
            }
        }
        //file << minGeneration<<"\n";
        generations--;
    }
    
    float val = Hill_Climbing(L, D, individual, funct);
    if(val < minFunc){
        file<<val<<"\n";
    }
    else{
        file<<minFunc<<"\n";
    }
    return minFunc;
    
}


int main()
{
    clock_t tStart = clock();
    srand(static_cast<unsigned int>(time(NULL)));
    //cout<<random01();
    int D = 5;
    int L = 16;
    int funct = 0;   //dixon Price
    //int funct = 1; // Rosenbrock
    //int funct = 2; // sphere
    //int funct = 3; // Rastrigin
    int generations = 1000;
    int pop_size = 100;
    
    ofstream f;
    f.open ("Gen_1000_Dixon_5.txt");
    
    for(int i = 0; i < 30; i++)
    {
        cout << genetic_algorithm(D, L, generations, funct, pop_size, f)<<endl;
        
    }
    
    /*for (int i = 0; i < population.size(); i++) {
        for (int j = 0; j < population[i].size(); j++)
            cout << population[i][j] << " ";
        cout << endl;
    }*/
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}
