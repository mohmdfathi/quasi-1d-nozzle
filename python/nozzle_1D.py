from params import Params
from grid   import grid_generation

if __name__ == "__main__":

    par = Params()

    x, A, dAdx = grid_generation( par )
    print("Grid size:", x.shape)
    print("Area min/max:", A.min(), A.max())

