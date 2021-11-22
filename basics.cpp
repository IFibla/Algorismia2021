#include <iostream>
#include <vector>
#include <set>

using Vertex = int;
using Adjacencia = set<Vertex>;
using Graf = vector<Adjacencia>;


bool esPIDS(const Graf& G, const set<Vertex>& D) {
    stack<int> S;

    S.push(0);

    while (not S.empty()) {
        int v = S.top();
        S.pop();

        
    }
}



int main() {
    cout << "Introdueix el Graf a analitzar i per finalitzar prem Ctrl+D" << endl;
    int v, a;
    while (cin >> v >> a) {
        Graf G(v);
        for (int i = 0; i < a ; ++i) {
            int x, y;
            cin >> x >> y;
            G[x].insert(y);
            G[y].insert(x);
        }
    }

    cout << "Introdueix el subconjunt D a comprovar i per finalitzar prem Ctrl+D" << endl;

    int v;
    set<Vertex> D();

    while (cin >> v) {
        D.insert(v);
    }

    if (esPIDS(G, D)) {
        if (esMPIDS(G, D)) 
            cout << "El conjunt es d'influencia positiva i tambe es minimal" << endl;
        else cout << "El conjunt es d'influecia positiva pero no es minimal" << endl;
    }

    else cout << "El conjunt no es d'influencia positiva" << endl;

}