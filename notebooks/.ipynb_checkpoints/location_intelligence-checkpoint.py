import osmnx as ox
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon
from scipy.interpolate import interp1d
import numpy as np
import networkx as nx
import geopandas as gpd

class LocationIntelligence:
    """ This calss provides a number of methods to work with osmnx apis and provide valueble output for location intelligence """
    def __init__(self, location, location_type, catchment_basis, basis_value):
        self.location = location
        self.location_type = location_type
        self.catchment_basis = catchment_basis
        self.basis_value = basis_value
        if self.location_type == 'address':
            self.network_graph = ox.graph_from_place(self.location, network_type='drive')
            self.network_graph = ox.project_graph(self.network_graph)
        elif self.location_type == 'point':
            if isinstance(center_point, tuple) == True:
                self.network_graph = ox.graph_from_point(center_point, dist=10000, dist_type='bbox', \
                                                         network_type='drive', simplify=True, retain_all=False, \
                                                         truncate_by_edge=False, clean_periphery=None, custom_filter=None)
                self.network_graph = ox.project_graph(self.network_graph)
            else:
                raise TypeError("Only tuple of (lat, lon) is allowed in case of point location_type")
        else:
            self.G = None
    
    def get_center_node(self):
        """ This method is to fetch the center node of the network based on the geocode result"""
        gdf_nodes = ox.graph_to_gdfs(self.network_graph, edges=False)
        x, y = gdf_nodes["geometry"].unary_union.centroid.xy
        center_node = ox.distance.nearest_nodes(self.network_graph, x[0], y[0])
        return center_node
        
    
    def create_catchment(self):
        """ Creates the catcment area around the address center node"""
        catchment_basis = self.catchment_basis
        basis_value = self.basis_value
        center_node = self.get_center_node()
        if catchment_basis == "time":
            if basis_value > 120:
                subgraph = None
            else:
                travel_speed = 20.0
                meters_per_minute = travel_speed * 1000 / 60  # km per hour to m per minute
                for _, _, _, data in self.network_graph.edges(data=True, keys=True):
                    data["time"] = data["length"] / meters_per_minute
                subgraph = nx.ego_graph(self.network_graph, center_node, radius=basis_value, distance="time")
        elif catchment_basis == "distance":
            if basis_value > 10:
                subgraph = None
            else:
                basis_value_m = basis_value * 1000.0
                subgraph = nx.ego_graph(self.network_graph, center_node, radius=basis_value_m, distance="length")
        if subgraph is not None:
            '''
            # uncomment if you want to plot the subgraph
            fig, ax = ox.plot_graph(
                subgraph, close=False, edge_color="red", edge_alpha=0.5, node_size=0, bgcolor='grey'
            )
            '''
            edge_buff=25 * basis_value
            node_buff=50 * basis_value
            node_points = [Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)]
            nodes_gdf = gpd.GeoDataFrame({"id": list(subgraph.nodes)}, geometry=node_points)
            nodes_gdf = nodes_gdf.set_index("id")
            edge_lines = []
            for n_fr, n_to in subgraph.edges():
                f = nodes_gdf.loc[n_fr].geometry
                t = nodes_gdf.loc[n_to].geometry
                edge_lookup = self.network_graph.get_edge_data(n_fr, n_to)[0].get("geometry", LineString([f, t]))
                edge_lines.append(edge_lookup)
            node_buffer = nodes_gdf.buffer(node_buff).geometry
            edge_buffer = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
            all_geom_buffer = list(node_buffer) + list(edge_buffer)
            node_edge_buffer = gpd.GeoSeries(all_geom_buffer).unary_union
            node_edge_buffer = Polygon(node_edge_buffer.exterior)
            simplified_polygon = node_edge_buffer.simplify(5)
            x = simplified_polygon.exterior.coords.xy[0]
            y = simplified_polygon.exterior.coords.xy[1]
            t = np.arange(len(x))
            ti = np.linspace(0, t.max(), 10 * t.size)
            xi = interp1d(t, x, kind='cubic')(ti)
            yi = interp1d(t, y, kind='cubic')(ti)
            final_polygon_shape = [[i,j] for i,j in zip(xi,yi)]
            final_polygon_shape = Polygon(final_polygon_shape)
            final_polygon_gdf = gpd.GeoDataFrame({"name":["catchment"]}, geometry=[final_polygon_shape])
            final_polygon_gdf.set_crs("epsg:{}".format(self.network_graph.graph['crs'].to_epsg()), inplace=True)
        else:
            final_polygon_gdf = None
        return final_polygon_gdf
        