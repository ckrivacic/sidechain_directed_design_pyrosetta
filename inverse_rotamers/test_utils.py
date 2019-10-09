from numeric import xyz_from_3d_array

def plot_3d(mobile,target,mobile_moved):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    fig = plt.Figure()
    ax = plt.axes(projection='3d')
    ax.scatter(*xyz_from_3d_array(mobile),c='blue')
    ax.scatter(*xyz_from_3d_array(target),c='green')
    ax.scatter(*xyz_from_3d_array(mobile_moved),c='red')
    plt.show()

