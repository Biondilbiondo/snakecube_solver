#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#define  I 0
#define _I 1
#define  J 2
#define _J 3
#define  K 4
#define _K 5

using std::cout;


class SnakeCube{
    private:
       //0 = stright
       //1 = angled
       int snake[64] = {0,0,1,1,0,1,1,1,0,0,1,1,0,1,1,0,
                        1,1,0,1,1,1,1,1,1,1,1,1,0,1,0,1,
                        1,1,1,1,1,0,1,0,0,1,1,1,1,0,0,1,
                        1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0};

       //Relative direction considering that the last cube is placed in the 
       //origin and oriented along z (k) axis. Each direction is
       //expressed as a versor oriented along or against the cartesian axis.
       int direction[64];

       //Absolute direction, first cube is considered to be oriented 
       //along z axis.
       int absolute_direction[64];

       //Absolute position, first cube is considered to be in the origin.
       int absolute_position[64][3];

       //Table that contains all possible results of a rotation of a versor
       int rot_table[6][6];

       //Just statistic for method solve(). Remove ASAP.
       int max_fold = 0;
       int solutions = 0;
       float explored_leafs = 0.0;

       //Just for testing performances
       double elapsed_rp = 0.0;
       double elapsed_rd = 0.0;
       double elapsed_gp = 0.0;
       double elapsed_gi = 0.0;
       

    public:
        SnakeCube( void ){
            for( int i = 0; i < 64; i++ )
                direction[ i ] = I;
            make_rotation_table();
            for( int i = 0; i < 6; i++ ){
                for( int j = 0; j < 6; j ++ )
                    printf("%d\t", rot_table[i][j] );
                printf("\n");
            }
            printf("\n");
            refresh_absolute_direction();
            refresh_absolute_position();
        }


        void make_rotation_table( void ){
            for( int v = 0; v < 6; v++ ){
                int vec[] = {0,0,0};
                if( v % 2 == 0 )
                    vec[ v/2 ] = 1;
                else
                    vec[ v/2 ] = -1;
                
                for( int w = 0; w < 6; w++ ){
                   int rot_matrix[3][3];
                   for( int j = 0; j < 3; j++ )
                       for( int k = 0; k < 3; k++ )
                            rot_matrix[j][k] = 0;
                   
                   switch( w ){
                       case I://+i rotation about y axis of pi/2
                           rot_matrix[0][2] = 1;
                           rot_matrix[1][1] = 1;
                           rot_matrix[2][0] = -1;
                           break;
                       case _I://-i rotation about y axis of -pi/2
                           rot_matrix[0][2] = -1;
                           rot_matrix[1][1] = 1;
                           rot_matrix[2][0] = 1;
                           break;
                       case J://+j rotation about x axis of pi/2
                           rot_matrix[0][0] = 1;
                           rot_matrix[1][2] = -1;
                           rot_matrix[2][1] = 1;
                           break;
                       case _J://-j
                           rot_matrix[0][0] = 1;
                           rot_matrix[1][2] = 1;
                           rot_matrix[2][1] = -1;
                           break;
                       case K://+k: identity
                           rot_matrix[0][0] = 1;
                           rot_matrix[1][1] = 1;
                           rot_matrix[2][2] = 1;
                           break;
                       case _K://-k rotation about x axis of pi
                           rot_matrix[0][0] = 1;
                           rot_matrix[1][1] = -1;
                           rot_matrix[2][2] = -1;
                           break;
                   }

                   int r[] = {0,0,0};
                   for( int a = 0; a < 3; a++ )
                       for( int b = 0; b < 3; b++ )
                            r[a] += rot_matrix[a][b] * vec[b];
                  
                   for( int a = 0; a < 3; a++ )
                       printf("%d\t", vec[a]);
                   printf("\n");

                   if( r[ 0 ] == 0 && r[ 2 ] == 0 ){
                       if( r[ 1 ] == 1 )
                           rot_table[v][w] = J;
                       else if( r[1] == -1 )
                           rot_table[v][w] = _J;
                   }
                   else if( r[ 1 ] == 0 && r[ 2 ] == 0 ){
                       if( r[ 0 ] == 1 )
                           rot_table[v][w] = I;
                       else if( r[0] == -1 )
                           rot_table[v][w] = _I;
                   }
                   else if( r[ 1 ] == 0 && r[ 0 ] == 0 ){
                       if( r[ 2 ] == 1 )
                           rot_table[v][w] = K;
                       else if( r[2] == -1 )
                           rot_table[v][w] = _K;
                   }
                   else rot_table[v][w] = 7;

                }
            }
        }

        int direction_id( int dir[] ){
            int res;
            //+i
            if( dir[0] == 1 && dir[1] == 0 && dir[2] == 0 )
                res = I;
            //+j
            if( dir[0] == 0 && dir[1] == 1 && dir[2] == 0 )
                res = J;
            //+k
            if( dir[0] == 0 && dir[1] == 0 && dir[2] == 1 )
                res = K;
            //-i
            if( dir[0] == -1 && dir[1] == 0 && dir[2] == 0 )
                res = _I;
            //-j
            if( dir[0] == 0 && dir[1] == -1 && dir[2] == 0 )
                res = _J;
            //-k
            if( dir[0] == 0 && dir[1] == 0 && dir[2] == -1 )
                res = _K;
            return res;
        }

        void refresh_absolute_direction( void ){
            double startTime = omp_get_wtime();
            int current_direction = K;
            for( int i = 0; i < 64; i++ ){
                if( get_element_type( i ) == 1 )
                   current_direction = rot_table[ direction[ i ] ][ current_direction ];
                absolute_direction[i] = current_direction; 
            }
        elapsed_rd +=  omp_get_wtime() - startTime ;
        }

        void refresh_absolute_position( void ){
            double startTime = omp_get_wtime();
            absolute_position[0][0] = 0;
            absolute_position[0][1] = 0;
            absolute_position[0][2] = 0;
            
            for( int i = 1; i < 64; i++ )
                for( int a = 0; a < 3; a++ ){
                    absolute_position[i][a] = absolute_position[i-1][a];
                    if( absolute_direction[ i ] / 2 == a ){
                        if( absolute_direction[ i ] % 2 == 0 )
                            absolute_position[i][a]++;
                        else
                            absolute_position[i][a]--;
                    }
                }
            
            elapsed_rp += omp_get_wtime() - startTime ;
        }

        int get_projection( int a, int topt = 64 ){
            double startTime = omp_get_wtime();
            int proj = 0;
            int min = 0;
            int max = 0;
                for( int i = 1; i < topt; i++ ){
                    if( absolute_position[i][a] < min )
                        min = absolute_position[i][a];
                    else if( absolute_position[i][a] > max )
                        max = absolute_position[i][a];
                }
                proj = max - min + 1;
            
            elapsed_gp +=  omp_get_wtime() - startTime;

            return proj;
        }
        
        bool no_intersection( int topt = 64 ){
            double startTime = omp_get_wtime();
            for( int j = 0; j < topt; j++ )
                for( int i = j + 1; i < topt; i++ )
                    if( absolute_position[i][0] == absolute_position[j][0] && \
                        absolute_position[i][1] == absolute_position[j][1] && \
                        absolute_position[i][2] == absolute_position[j][2] )
                        return 0;
            elapsed_gi +=  omp_get_wtime() - startTime;
            return 1;
        }

        bool switch_element( int el, int rot_id ){
            if( rot_id < 0 || rot_id > 4 )
                return 1;
            if( get_element_type( el ) == 0 )
                return 1;
            
            direction[ el ] = rot_id;
            
            refresh_absolute_direction();
            refresh_absolute_position();
            return 0;
        }

        void print_absolute_position( void ){
            for( int j = 0; j < 8; j++ ){
                for( int i = j*8; i < (j+1)*8; i++ )
                    cout << std::setw(10) << i << "\t";
                cout << "\n";
                for( int i = j*8; i < (j+1)*8; i++ ){
                    cout << '(' << std::setw(2) << absolute_position[i][0]; 
                    cout << ';' << std::setw(2) << absolute_position[i][1]; 
                    cout << ';' << std::setw(2) << absolute_position[i][2] << ')'; 
                    cout << "\t";
                }
                cout << "\n\n";
            }
        }

        void print_absolute_direction( void ){
            for( int j = 0; j < 8; j++ ){
                for( int i = j*8; i < (j+1)*8; i++ )
                    cout << std::setw(2) << i << "\t";
                cout << "\n";
                for( int i = j*8; i < (j+1)*8; i++ ){
                    switch( absolute_direction[ i ] ){
                        case 0:
                            cout << "+i";
                            break;
                        case 1:
                            cout << "-i";
                            break;
                        case 2:
                            cout << "+j";
                            break;
                        case 3:
                            cout << "-j";
                            break;
                        case 4:
                            cout << "+k";
                            break;
                        case 5:
                            cout << "-k";
                            break;
                        default:
                            cout << "NK";
                            break;
                    }
                    cout << "\t";
                }
                cout << "\n\n";
            }
        }

        int get_element_type( int i ){
            if( i >= 64 || i < 0 ){
                return -1;
            }
            return snake[ i ];
        }

        bool solve( int start = 0 ){
            if( get_projection( 0, start ) > 4 || \
                get_projection( 1, start ) > 4 || \
                get_projection( 2, start ) > 4 ){
                int cnt = 0;
                for( int j = start + 1; j < 64; j++ )
                    cnt += get_element_type( j );
                explored_leafs += pow( 4.0, cnt );
                return 1;
            }
            if( no_intersection(start) == 0 ){
                int cnt = 0;
                for( int j = start + 1; j < 64; j++ )
                    cnt += get_element_type( j );
                explored_leafs += pow( 4.0, cnt );
                return 1;
            }
            
            int next_el = -1;
            for( int i = start + 1; i < 64; i++ )
                if( get_element_type( i ) == 1 ){
                    next_el = i;
                    break;
                }
            
            if( next_el == -1 ){
                cout << "New maximum depth reached in folding: " << max_fold << '\n';
                cout << "Explored leafs " << 100.00 * explored_leafs/1.24e27 << "\n\n";
                
                if( get_projection( 0 ) > 4 || \
                    get_projection( 1 ) > 4 || \
                    get_projection( 2 ) > 4 ){
                    int cnt = 0;
                    for( int j = start + 1; j < 64; j++ )
                        cnt += get_element_type( j );
                    explored_leafs += pow( 4.0, cnt );
                    return 1;
                }
                if( no_intersection() == 0 ){
                    int cnt = 0;
                    for( int j = start + 1; j < 64; j++ )
                        cnt += get_element_type( j );
                    explored_leafs += pow( 4.0, cnt );
                    return 1;
                }
                
                solutions++;
                cout << "Found a solution!!!\n";
                return 1;
                return 0;
            }

            if( next_el > max_fold ){
                max_fold = next_el;
                cout << "New maximum depth reached in folding: " << max_fold << '\n';
                cout << "Explored leafs " << 100.00 * explored_leafs/1.24e27 << "\n\n";
            }
            
            for( int a = 0; a < 4; a++ ){
                switch_element( next_el, a );
                if( solve( next_el ) == 0 )
                    return 0;
            }

            if( start == 0 ){
                cout << "Total solutions: " << solutions << '\n';
            }
            return 1;
        }
};

int main( void ){
    SnakeCube sc;
    sc.print_absolute_direction();
    sc.print_absolute_position();
    for( int i = 0; i < 3; i++ )
        cout << "Projection along " << i << ": " << sc.get_projection(i) << "\n";
    if( sc.no_intersection() == 0 )
        cout << "An intersection is present.";
    int cnt = 0;
    for( int i = 0; i < 64; i++ )
        cnt += sc.get_element_type( i );
    cout << "Number of angles: " << cnt << "\n";

    if( sc.solve() == 1 )
        cout << "Folding failed.\n";
    else{
        cout << "Folding finished\n";
        sc.print_absolute_direction();
        sc.print_absolute_position();
        for( int i = 0; i < 3; i++ )
            cout << "Projection along " << i << ": " << sc.get_projection(i) << "\n";
        if( sc.no_intersection() == 0 )
            cout << "An intersection is present.";
        cnt = 0;
        for( int i = 0; i < 64; i++ )
            cnt += sc.get_element_type( i );
        cout << "Number of angles: " << cnt << "\n";
    }
}
