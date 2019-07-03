//
//  main.cpp
//  WOW
//
//  Created by yacymelo on 2019/7/1.
//  Copyright © 2019 yacymelo. All rights reserved.
//

#include <iostream>
#include <climits>
#include "include/Graph.h"

using namespace std;

typedef unsigned int depth_type;

depth_type MAX_DIST = std::numeric_limits<depth_type>::max();

class Vetex_Prop {
public:
    depth_type depth;
    //int parent;
    int id;
    bool active;
public:
    Vetex_Prop() {
        depth = MAX_DIST;
        //parent = -1;
        id = -1;
        active = false;
    }
    bool operator != (const Vetex_Prop& p) {
        return (this->depth != p.depth);
    }
    
    void operator = (const Vetex_Prop& p) {
        depth = p.depth;
        //parent = p.parent;
        id = p.id;
        active = p.active;
    }
    
    friend std::ostream &operator<<(std::ostream &outstream, const Vetex_Prop & val)
    {
        outstream << val.depth;
        return outstream;
    }
    
    void Process_Edge (Vetex_Prop& Res) {
        Res.depth = depth+1;
    }
    
    void Reduce (Vetex_Prop& Res) {
        if (depth > Res.depth)
            depth = Res.depth;
    }
    
    void Apply(Vetex_Prop& tmp)  {
        if (depth > tmp.depth){
            depth = tmp.depth;
            active = true;
        }
    }
    
    void Reset()  {
        active = false;
    }
    
    void ResetTmp()  {
        depth = MAX_DIST;
        active = false;
    }
    
    void RootInti () {
        depth = 0;
        active =true;
    }
//    void do_every_iteration(int iteration_number) {
//        current_depth++;
//    }
};

void run_bfs(char * filename, int v) {
    Graph<Vetex_Prop> G;
    G.LoadGraph(filename);
    int root = v;
    G.SetRoot(root);
    G.GRun_Syn_Pro(UNTIL_CONVERGENCE);
}

int main(int argc, char * argv[]) {
    if (argc < 3) {
        printf("Correct format: %s A.mtx source_vertex (1-based index)\n", argv[0]);
        return 0;
    }
    int source_vertex = atoi(argv[2]);
    run_bfs(argv[1], source_vertex);
    return 0;
}

