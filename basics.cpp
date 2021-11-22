#include <iostream>
#include <vector>
#include <set>
#include <stack>
using namespace std;

using Adjacencia = set<int>;
using Graf = vector<Adjacencia>;


bool isPositiveInfluenceDominatingSet(const Graf& G, const set<int>& D, vector<bool>& visitat) {

    stack<int> S;

    S.push(0);
    bool esPIDS = true;                             // Comencem assumint que ho es

    while (not S.empty() and esPIDS) {

        int v = S.top();
        S.pop();

        int nAdjacenciesD = 0;

        if (not visitat[v]) {

            visitat[v] = true;

            for (int aresta : G[v]) {

                std::set<int>::iterator it = D.find(aresta);
                if (it != D.end()) ++nAdjacenciesD;
                S.push(aresta);

            }

            int nAdjacenciesTotals = G[v].size();

            double proporcioAdjaD = double(nAdjacenciesD)/double(nAdjacenciesD);
            cout << "Proporcio graus adj vertex " << v << " " << proporcioAdjaD << " " << nAdjacenciesD << endl;
            if (proporcioAdjaD < 0.5) esPIDS = false;

        }
    }
    return esPIDS;
}

int main() {

    cout << "Lectura del Graf G = (V,E)" << endl;
    cout << "Insereix el nombre de vertexs n = |V| i el nombre d'arestes m = |E|" << endl;

    int n, m;
    cin >> n >> m;

    Graf G(n);
    cout << "Introdueix m arestes de la forma v u on v,u son vertexs de V" << endl;

    int v, u;
    while (m--) {
        cin >> v >> u;
        G[v].insert(u);
        G[u].insert(v);
    }

    cout << "Introdueix el nombre de vertexs del conjunt D i seguidament els vertexs que hi pertanyen" << endl;

    int nD;
    cin >> nD;
    set<int> D;
    
    while (nD--) {
        cin >> v;
        D.insert(v);
    }

    vector<bool> visitat(n, false);
/*
    for (int i = 0; i < n; ++i) {
        cout << "VERTEX " << i << endl;
        cout << "Arestes ";
        for (int aresta : G[i]) {
            cout << aresta << " ";
        }
        cout << endl;
    }
*/

    cout << "Comprovant si D es PIDS..." << endl;

    if (isPositiveInfluenceDominatingSet(G, D, visitat)) {
        //if (esMPIDS(G, D)) 
            cout << "El conjunt es d'influencia positiva" << endl;
        //else cout << "El conjunt es d'influecia positiva pero no es minimal" << endl;
    }

    else cout << "El conjunt no es d'influencia positiva" << endl;

}