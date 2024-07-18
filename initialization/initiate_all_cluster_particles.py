import initialize_isotropic_plummers as IIP
import sys 
import tstrippy



if __name__=="__main__":
    cluster_num     =   int(sys.argv[1])
    NP              =   int(sys.argv[2])

    cluster_names = tstrippy.Parsers.baumgardtMWGCs().data['Cluster'][:]
    N_clusters = len(cluster_names)

    assert cluster_num < N_clusters, "Cluster number out of range"

    GCname = cluster_names[cluster_num]
    IIP.main(GCname, NP)