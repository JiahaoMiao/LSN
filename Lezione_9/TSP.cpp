#include "TSP.h"

//Constructor of the path class
path::path(std::vector<int> path):_path(std::move(path)){
    if(!CheckPath(_path)){
        std::cerr << "error: the path is not valid\n";
        exit(-1);
    }
}

//check if all IDs are unique
bool path::CheckPath(const std::vector<int>& vec){
    std::unordered_set<int> labels;
    
    // Controlla se il vettore contiene 33 città, la prima città è 0 per tutti
    // if (vec.size() != 33) {
    //     std::cerr << "error: the path does not contain 34 cities\n";
    //     return false;
    // }
    
    // Controlla se il numero è già presente nel set
    if (!std::all_of(vec.begin(), vec.end(), [&](int num){ return labels.insert(num).second; }) ) {
        std::cerr << "error: the path contains duplicate cities\n";
        return false;
    }
    
    return true;
}


//Calculate the distance between two cities using the L1 norm
inline double path::distanceL1(city c, city d){
    return fabs(c.getX()-d.getX())+fabs(c.getY()-d.getY());
}

//Calculate the distance between two cities using the L2 norm
inline double path::distanceL2(city c, city d){
    return sqrt(pow(c.getX()-d.getX(),2)+pow(c.getY()-d.getY(),2));
}

double path::LengthL1(const std::vector<city>& cities){
    if (_path.empty()) return 0.0;
    //initialize the length with the distance between the first and the last city and the first and second city
    //the first city is the one with ID 0 and the last one is the one with ID _path[_path.size()-1] 
    double length = {distanceL1(cities[0],cities[_path[_path.size()-1]])+distanceL1(cities[0],cities[_path[0]])};
    //add the distance between the current city and the next one
    for(unsigned int i{}; i < _path.size()-1; i++){
        length += distanceL1(cities[_path[i]],cities[_path[i+1]]);
    }
    return length;
}

double path::LengthL2(const std::vector<city>& cities){
    if (_path.empty()) return 0.0;
    //initialize the length with the distance between the first and the last city and the first and second city
    double length = {distanceL2(cities[0],cities[_path[_path.size()-1]]) + distanceL2(cities[0],cities[_path[0]])};
    for(unsigned int i{}; i < _path.size()-1; i++){
        length += distanceL2(cities[_path[i]],cities[_path[i+1]]);
    }
    return length;
}

//Print the path
void path::PrintPath(){
    std::cout << "[0, ";
    std::for_each(_path.begin(), _path.end(), [](int i) { std::cout << i << ", "; });
    std::cout << "\b\b]\n";
}

//shift of m contiguous cities of n positions starting from the i-th city
void path::shift(int i,int m, int n){
    // apply swap m times
    // i--; //convention: because i work with a vector with 0 as first element and i work with the rest of the cities
    for(int j = 0; j < m; j++){
        //swap the i+j-th city with the (i+j+n)-th city
        std::swap(_path[(i+j)%_path.size()],_path[(i+j+n)%_path.size()]);
    }
}

//Permutation of m cities with other different m cities
//start from the "start"-th city and separate the two groups of cities by "sep" cities
void path::Permutation(int m,int start,int sep){
    m %= _path.size()/2; //m must be less than the half of the path
    int size = _path.size();
    int diff = size - m;
    sep %= diff; //separation must be less than the difference between the length and m, so that i have m cities to swap with 
    for(int j = 0; j < m; j++){
        //swap the (start+j)-th city with the (start+m+sep+j)-th city
        std::swap(_path[(start+j)%_path.size()],_path[(start+m+sep+j)%_path.size()]);
    }
}

//inversion of the order of m contiguous cities starting from the i-th city
void path::inversion(int m,int i){
    m %= _path.size(); //m must be less than the length of the path
    int size = _path.size();
    for(int j = 0; j < m/2; j++){
        //swap the (i+j)-th city with the (i+m-j-1)-th city
        // i, i+1, i+2, ..., i+m-1 (are m cities)
        std::swap(_path[(i+j)%size],_path[(i+m-1-j)%size]);
    }
}

//Constructor       
TSP::TSP(){
    _rnd.initialize();
}

TSP::TSP(std::vector<city> cities):_cities(std::move(cities)){
    _rnd.initialize();
}

TSP::TSP(std::vector<city> cities,std::vector<path> paths):_cities(std::move(cities)),_paths(std::move(paths)){
    _rnd.initialize();
}

//pair permutation of cities
path& TSP::PairPermutation(path& way){
    int i = _rnd.Rannyu(0,_cities.size()-1);
    int j = _rnd.Rannyu(0,_cities.size()-1);
    while(i==j){
        j = _rnd.Rannyu(0,_cities.size()-1);
    }
    way.swap(i,j);
    return way;
}

//shift of m contiguous cities of n positions starting from the i-th city
path& TSP::Shift(path& way,int m, int n) {
    way.shift(_rnd.Rannyu(0,_cities.size()-1), m, n);
    return way;
}

//permutation among m contiguous cities with other m contiguous cities
path& TSP::Permutation(path& way,int m){
    int start = _rnd.Rannyu(0,_cities.size()-1);
    int sep = _rnd.Rannyu(0,_cities.size()-1);
    way.Permutation(m,start,sep);
    return way;
}

//inversion of the order of m contiguous cities starting from the i-th city
path& TSP::Inversion(path& way,int m){
    int start = _rnd.Rannyu(0,_cities.size()-1);
    way.inversion(m,start);
    return way;
}

//Selection
int TSP::Selection(){
    //use M*rand^p to select the path
    //p > 1 to have greater probability to select the first half of the paths(if the paths are sorted in ascending order) 
    //p < 1 to have greater probability to select the second half of the paths(if the paths are sorted in descending order)
    //p = 2 is equivalent of having 70.7% of the probability to select the first half of the paths 
    //p = 3 is equivalent of having 79.3% of the probability to select the first half of the paths
    //p = 4 is equivalent of having 84.1% of the probability to select the first half of the paths
    return _paths.size()*pow(_rnd.Rannyu(),3);
}

//New generation of paths
void TSP::NewGeneration(){
    std::vector<path> new_gen(_paths.size());
    //CrossOver
    
    //using L1 norm to sort the paths in ascending order (shortest path first)
    // sort(_paths.begin(),_paths.end(),[&](path& a, path& b){return a.LengthL1(_cities) < b.LengthL1(_cities);});
    //using L2 norm to sort the paths in ascending order (shortest path first)
    sort(_paths.begin(),_paths.end(),[&](path& a, path& b){return a.LengthL2(_cities) < b.LengthL2(_cities);});
   
    unsigned int start{0u};
    //Elitism: if _paths.size() is even, only the best path is copied to the new generation
    //if _paths.size() is odd, the best two paths are copied to the new generation
    if(_paths.size()%2!=0){
        new_gen[0] = _paths[0];
        start++; //start from the second element
    }else{
        new_gen[0] = _paths[0];
        new_gen[1] = _paths[1];
        start+=2u; //start from the third element
    }

    int cut{},i{},j{};
    int N = _cities.size();
    path way;
    for(unsigned int k{start};k<_paths.size();k+=2){
        
        //choose random path
        i = Selection(); //first parent
        j = Selection(); //second parent
        while(i==j){
            j = Selection(); //ensure that the two parents are different
        }
        std::vector<int> child1(N-1);//the first city is always 0, so i omit it
        std::vector<int> child2(N-1);
        
        //Crossover
        if(_rnd.Rannyu()<0.60){//60% of the probability to perform crossover
            //choose random cut point
            cut = _rnd.Rannyu(0,N-1);
            
            // copy elements before cut from the first parent to the child1
            std::copy_n(_paths[i].getPath().begin(),cut,child1.begin());
            std::copy_n(_paths[j].getPath().begin(),cut,child2.begin());

            // Add the remaining elements of the second parent to the child if they are not already present
            // same for the seccond child
            for(int l{}, m1{cut}, m2{cut}; l<N-1; l++){
                //if the value is not already present in the child, add it
                if(std::count(child1.begin(),child1.end(),_paths[j].getPath()[l])==0){
                    child1[m1++] = _paths[j].getPath()[l];
                }
                if(std::count(child2.begin(),child2.end(),_paths[i].getPath()[l])==0){
                    child2[m2++] = _paths[i].getPath()[l];
                }
            }

            new_gen[k] = path(std::move(child1));
            new_gen[k+1] = path(std::move(child2));
        }else{
            //if the crossover is not performed, copy the parents to the new generation
            new_gen[k] = _paths[i];
            new_gen[k+1] = _paths[j];
        }
    }
    //Mutations
    for(auto& way: new_gen){
        Mutate(way);
    }
    
    // Assign the new generation of paths
    _paths = new_gen;

    //Sort the paths in ascending order (shortest path first)    
    // sort(_paths.begin(),_paths.end(),[&](path& a, path& b){return a.LengthL1(_cities) < b.LengthL1(_cities);});
    sort(_paths.begin(),_paths.end(),[&](path& a, path& b){return a.LengthL2(_cities) < b.LengthL2(_cities);});
}

//Mutations
path& TSP::Mutate(path& way){
    //10% of the probability to perform mutation
    if(_rnd.Rannyu()< 0.2){
        // std::cout << "swap happened\n";
        way.swap(_rnd.Rannyu(0,_cities.size()-1),_rnd.Rannyu(0,_cities.size()-1));
        // if(!way.CheckPath(way.getPath())){
        //     std::cerr << "error: the path is not valid\n";
        //     exit(-1);
        // }
    }
    if(_rnd.Rannyu()<0.1){
        // std::cout << "shift happened\n";
        way.shift(_rnd.Rannyu(0,_cities.size()-1),_rnd.Rannyu(1,_cities.size()-1),_rnd.Rannyu(0,_cities.size()-1));
        // if(!way.CheckPath(way.getPath())){
        //     std::cerr << "error: the path is not valid\n";
        //     exit(-1);
        // }
    }
    if(_rnd.Rannyu()<0.1){
        // std::cout << "permutation happened\n";
        // m>1 to avoid the permutation of a single city(swap)
        way.Permutation(_rnd.Rannyu(1,_cities.size()-1),_rnd.Rannyu(0,_cities.size()-1),_rnd.Rannyu(0,_cities.size()-1));
        // if(!way.CheckPath(way.getPath())){
        //     std::cerr << "error: the path is not valid\n";
        //     exit(-1);
        // }
    }
    if(_rnd.Rannyu()<0.15){
        way.inversion(_rnd.Rannyu(0,_cities.size()-1),_rnd.Rannyu(0,_cities.size()-1)); 
        // if(!way.CheckPath(way.getPath())){
        //     std::cerr << "error: the path is not valid\n";
        //     exit(-1);
        // }
    }
    return way;
}   

//Print the path with length
void TSP::Print(){
    for(auto& p: _paths){
        p.PrintPath();
        std::cout << "Length L1: " << p.LengthL1(_cities) << "  L2: " << p.LengthL2(_cities) << "\n";
    }
}

//Print the best path
void TSP::PrintBest(){
    _paths[0].PrintPath();
    std::cout << "Length L1: " << _paths[0].LengthL1(_cities) << "  L2: " << _paths[0].LengthL2(_cities) << "\n";
}

//Average of the best half of the paths
double TSP::AverageLength(){
    double sum{};
    for(unsigned int i{}; i < _paths.size()/2; i++){
        sum += _paths[i].LengthL2(_cities);
    }
    return sum/double(_paths.size()/2);
}
//return the best length of the current population
double TSP::BestLength(){
    return _paths[0].LengthL2(_cities);
}

//Print the best path to file
void TSP::BestPath(const std::string& filename){
    std::ofstream out(filename);
    if(!out){
        std::cerr << "error: unable to open the file" << filename << "\n";
        exit(-1);
    }
    //Position of the first city
    out <<_cities[0].getX() << " " << _cities[0].getY() << "\n";
    // _paths[0] is the best path (already sorted in ascending order)
    for(auto i: _paths[0].getPath()){
        out << _cities[i].getX() << " " << _cities[i].getY() << "\n";
    }
    //Position of the first city to close the path
    out << _cities[0].getX() << " " << _cities[0].getY() << "\n";
    
    out.close();
}
