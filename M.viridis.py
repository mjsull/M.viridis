#!/usr/bin/env python
# draw_coloured_tree.py   Written by: Mitchell Sullivan   mjsull@gmail.com
# organisation: Icahn School of Medicine - Mount Sinai
# Version 0.0.1 2016.01.19
# License: GPLv3

from ete3 import Tree, RectFace, AttrFace, TextFace
import sys
from ete3 import NodeStyle
from ete3 import TreeStyle
import argparse
import datetime
import subprocess
import os
import struct

def colorstr(rgb):
    return "#%02x%02x%02x" % (rgb[0],rgb[1],rgb[2])

def strtorgb(rgbstr):
    return struct.unpack('BBB',rgbstr[1:].decode('hex'))

# take a hue, saturation and lightness value and return a RGB hex string
def hsl_to_str(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
    m = l - c/2
    if h < 60:
        r, g, b = c + m, x + m, 0 + m
    elif h < 120:
        r, g, b = x + m, c+ m, 0 + m
    elif h < 180:
        r, g, b = 0 + m, c + m, x + m
    elif h < 240:
        r, g, b, = 0 + m, x + m, c + m
    elif h < 300:
        r, g, b, = x + m, 0 + m, c + m
    else:
        r, g, b, = c + m, 0 + m, x + m
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return "#%02x%02x%02x" % (r,g,b)

# main function for drawing tree
def draw_tree(the_tree, colour, back_color, label, out_file, the_scale, extend, bootstrap, group_file, grid_options, the_table, pres_abs, circular):
    t = Tree(the_tree, quoted_node_names=True)
#    t.ladderize()
    font_size = 8
    font_type = 'Heveltica'
    font_gap = 3
    font_buffer = 10
    o = t.get_midpoint_outgroup()
    t.set_outgroup(o)
    the_leaves = []
    for leaves in t.iter_leaves():
        the_leaves.append(leaves)
    groups = {}
    num = 0
    # set cutoff value for clades as 1/20th of the distance between the furthest two branches
    # assign nodes to groups
    last_node = None
    ca_list = []
    if not group_file is None:
        style = NodeStyle()
        style['size'] = 0
        style["vt_line_color"] = '#000000'
        style["hz_line_color"] = '#000000'
        style["vt_line_width"] = 1
        style["hz_line_width"] = 1
        for n in t.traverse():
            n.set_style(style)
        with open(group_file) as f:
            group_dict = {}
            for line in f:
                group_dict[line.split()[0]] = line.split()[1]
        for node in the_leaves:
            i = node.name
            for j in group_dict:
                if j in i:
                    if group_dict[j] in groups:
                        groups[group_dict[j]].append(i)
                    else:
                        groups[group_dict[j]] = [i]
        coloured_nodes = []
        for i in groups:
            the_col = i
            style = NodeStyle()
            style['size'] = 0
            style["vt_line_color"] = the_col
            style["hz_line_color"] = the_col
            style["vt_line_width"] = 2
            style["hz_line_width"] = 2
            if len(groups[i]) == 1:
                ca = t.search_nodes(name=groups[i][0])[0]
                ca.set_style(style)
                coloured_nodes.append(ca)
            else:
                ca = t.get_common_ancestor(groups[i])
                ca.set_style(style)
                coloured_nodes.append(ca)
                tocolor = []
                for j in ca.children:
                    tocolor.append(j)
                while len(tocolor) > 0:
                    x = tocolor.pop(0)
                    coloured_nodes.append(x)
                    x.set_style(style)
                    for j in x.children:
                        tocolor.append(j)
            ca_list.append((ca, the_col))
        if back_color:
            # for each common ancestor node get it's closest common ancestor neighbour and find the common ancestor of those two nodes
            # colour the common ancestor then add it to the group - continue until only the root node is left
            while len(ca_list) > 1:
                distance = float('inf')
                for i, col1 in ca_list:
                    for j, col2 in ca_list:
                        if not i is j:
                            parent = t.get_common_ancestor(i, j)
                            getit = True
                            the_dist = t.get_distance(i, j)
                            if the_dist <= distance:
                                distance = the_dist
                                the_i = i
                                the_j = j
                                the_i_col = col1
                                the_j_col = col2
                ca_list.remove((the_i, the_i_col))
                ca_list.remove((the_j, the_j_col))
                rgb1 = strtorgb(the_i_col)
                rgb2 = strtorgb(the_j_col)
                rgb3 = ((rgb1[0] + rgb2[0])/2, (rgb1[1] + rgb2[1])/2, (rgb1[2] + rgb2[2])/2)
                new_col = colorstr(rgb3)
                new_node = t.get_common_ancestor(the_i, the_j)
                the_col = new_col
                style = NodeStyle()
                style['size'] = 0
                style["vt_line_color"] = the_col
                style["hz_line_color"] = the_col
                style["vt_line_width"] = 2
                style["hz_line_width"] = 2
                new_node.set_style(style)
                coloured_nodes.append(new_node)
                ca_list.append((new_node, new_col))
                for j in new_node.children:
                    tocolor.append(j)
                while len(tocolor) > 0:
                    x = tocolor.pop(0)
                    if not x in coloured_nodes:
                        coloured_nodes.append(x)
                        x.set_style(style)
                        for j in x.children:
                            tocolor.append(j)
    elif colour:
        distances = []
        for node1 in the_leaves:
            for node2 in the_leaves:
                if node1 != node2:
                    distances.append(t.get_distance(node1, node2))
        distances.sort()
        clade_cutoff = distances[len(distances)/4]
        for node in the_leaves:
            i = node.name
            if not last_node is None:
                if t.get_distance(node, last_node) <= clade_cutoff:
                    groups[group_num].append(i)
                else:
                    groups[num] = [num, i]
                    group_num = num
                    num += 1
            else:
                groups[num] = [num, i]
                group_num = num
                num += 1
            last_node = node
        for i in groups:
            num = groups[i][0]
            h = num * 360/len(groups)
            the_col = hsl_to_str(h, 0.5, 0.5)
            style = NodeStyle()
            style['size'] = 0
            style["vt_line_color"] = the_col
            style["hz_line_color"] = the_col
            style["vt_line_width"] = 2
            style["hz_line_width"] = 2
            if len(groups[i]) == 2:
                ca = t.search_nodes(name=groups[i][1])[0]
                ca.set_style(style)
            else:
                ca = t.get_common_ancestor(groups[i][1:])
                ca.set_style(style)
                tocolor = []
                for j in ca.children:
                    tocolor.append(j)
                while len(tocolor) > 0:
                    x = tocolor.pop(0)
                    x.set_style(style)
                    for j in x.children:
                        tocolor.append(j)
            ca_list.append((ca, h))
        # for each common ancestor node get it's closest common ancestor neighbour and find the common ancestor of those two nodes
        # colour the common ancestor then add it to the group - continue until only the root node is left
        while len(ca_list) > 1:
            distance = float('inf')
            got_one = False
            for i, col1 in ca_list:
                for j, col2 in ca_list:
                    if not i is j:
                        parent = t.get_common_ancestor(i, j)
                        getit = True
                        for children in parent.children:
                            if children != i and children != j:
                                getit = False
                                break
                        if getit:
                            the_dist = t.get_distance(i, j)
                            if the_dist <= distance:
                                distance = the_dist
                                the_i = i
                                the_j = j
                                the_i_col = col1
                                the_j_col = col2
                                got_one = True
            if not got_one:
                break
            ca_list.remove((the_i, the_i_col))
            ca_list.remove((the_j, the_j_col))
            new_col = (the_i_col + the_j_col) / 2
            new_node = t.get_common_ancestor(the_i, the_j)
            the_col = hsl_to_str(new_col, 0.5, 0.3)
            style = NodeStyle()
            style['size'] = 0
            style["vt_line_color"] = the_col
            style["hz_line_color"] = the_col
            style["vt_line_width"] = 2
            style["hz_line_width"] = 2
            new_node.set_style(style)
            ca_list.append((new_node, new_col))
    # if you just want a black tree
    else:
        style = NodeStyle()
        style['size'] = 0
        style["vt_line_color"] = '#000000'
        style["hz_line_color"] = '#000000'
        style["vt_line_width"] = 1
        style["hz_line_width"] = 1
        for n in t.traverse():
            n.set_style(style)
    color_list = [(240,163,255),(0,117,220),(153,63,0),(76,0,92),(25,25,25),(0,92,49),(43,206,72),(255,204,153),
                  (128,128,128),(148,255,181),(143,124,0),(157,204,0),(194,0,136),(0,51,128),(255,164,5),(255,168,187),
                  (66,102,0),(255,0,16),(94,241,242),(0,153,143),(224,255,102),(116,10,255),(153,0,0),(255,255,128),
                  (255,255,0),(255,80,5), (0, 0, 0), (50, 50, 50)]
    up_to_colour = {}
    ts = TreeStyle()
    column_list = []
    width_dict = {}
    if not grid_options is None:
        colour_dict = {}
        type_dict = {}
        min_val_dict = {}
        max_val_dict = {}
        leaf_name_dict = {}
        header_count = 0
        the_columns = {}
        if grid_options == 'auto':
            with open(the_table) as f:
                headers = f.readline().rstrip().split('\t')[1:]
                for i in headers:
                    the_columns[i] = [i]
                    type_dict[i] = 'colour'
                    colour_dict[i] = {'empty':'#FFFFFF'}
                    width_dict[i] = 20
                    up_to_colour[i] = 0
                    column_list.append(i)
        else:
            with open(grid_options) as g:
                for line in g:
                    if line.startswith('H'):
                        name, type, width = line.rstrip().split('\t')[1:]
                        if name in the_columns:
                            the_columns[name].append(name + '_' + str(header_count))
                        else:
                            the_columns[name] = [name + '_' + str(header_count)]
                        width = int(width)
                        name = name + '_' + str(header_count)
                        header_count += 1
                        colour_dict[name] = {'empty':'#FFFFFF'}
                        type_dict[name] = type
                        width_dict[name] = width
                        column_list.append(name)
                        up_to_colour[name] = 0
                        min_val_dict[name] = float('inf')
                        max_val_dict[name] = 0
                    elif line.startswith('C'):
                        c_name, c_col = line.rstrip().split('\t')[1:]
                        if not c_col.startswith('#'):
                            c_col = colorstr(map(int, c_col.split(',')))
                        colour_dict[name][c_name] = c_col
        val_dict = {}
        with open(the_table) as f:
            headers = f.readline().rstrip().split('\t')[1:]
            column_no = {}
            for num, i in enumerate(headers):
                if i in the_columns:
                    column_no[num] = i
            for line in f:
                name = line.split('\t')[0]
                leaf_name = None
                for n in t.traverse():
                    if n.is_leaf():
                        if name.split('.')[0] in n.name:
                            leaf_name = n.name
                if leaf_name is None:
                    continue
                else:
                    leaf_name_dict[leaf_name] = name
                vals = line.rstrip().split('\t')[1:]
                if name in val_dict:
                    sys.exit('Duplicate entry found in table.')
                else:
                    val_dict[name] = {}
                for num, val in enumerate(vals):
                    if num in column_no and val != '':
                        for q in the_columns[column_no[num]]:
                            column_name = q
                            if type_dict[column_name] == 'colour':
                                val_dict[name][column_name] = val
                                if not val in colour_dict[column_name]:
                                    colour_dict[column_name][val] = colorstr(color_list[up_to_colour[column_name] % len(color_list)])
                                    up_to_colour[column_name] += 1
                            elif type_dict[column_name] == 'text':
                                val_dict[name][column_name] = val
                            elif type_dict[column_name] == 'colour_scale_date':
                                year, month, day = val.split('-')
                                year, month, day = int(year), int(month), int(day)
                                the_val = datetime.datetime(year, month, day, 0, 0, 0) - datetime.datetime(1970, 1, 1, 0, 0, 0)
                                val_dict[name][column_name] = the_val.total_seconds()
                                if the_val.total_seconds() < min_val_dict[column_name]:
                                    min_val_dict[column_name] = the_val.total_seconds()
                                if the_val.total_seconds() > max_val_dict[column_name]:
                                    max_val_dict[column_name] = the_val.total_seconds()
                            elif type_dict[column_name] == 'colour_scale':
                                the_val = float(val)
                                val_dict[name][column_name] = the_val
                                if the_val < min_val_dict[column_name]:
                                    min_val_dict[column_name] = the_val
                                if the_val > max_val_dict[column_name]:
                                    max_val_dict[column_name] = the_val
                            else:
                                sys.exit('Unknown column type')
        if not out_file is None:
            new_desc = open(out_file + '.new_desc', 'w')
        else:
            new_desc = open('viridis.new_desc', 'w')
        ts.legend_position = 3
        leg_column = 0
        for num, i in enumerate(column_list):
            nameF = TextFace(font_gap * ' ' + i.rsplit('_',1)[0] + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True)
            nameF.rotation = -90
            ts.aligned_header.add_face(nameF, column=num+1)
            new_desc.write('H\t' + i.rsplit('_', 1)[0] + '\t' + type_dict[i] + '\t' + str(width_dict[i]) + '\n')
            x = num * 200
            if type_dict[i] == 'colour':
                ts.legend.add_face(TextFace(font_gap * ' ' + i.rsplit('_',1)[0] + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                ts.legend.add_face(RectFace(width_dict[i], 20, '#FFFFFF', '#FFFFFF'), column=leg_column)
                for num2, j in enumerate(colour_dict[i]):
                    new_desc.write('C\t' + j + '\t' + colour_dict[i][j] + '\n')
                    ts.legend.add_face(TextFace(font_gap * ' ' + j + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                    ts.legend.add_face(RectFace(width_dict[i], 20, colour_dict[i][j], colour_dict[i][j]), column=leg_column)
                leg_column += 2
            elif type_dict[i] == 'colour_scale':
                ts.legend.add_face(TextFace(font_gap * ' ' + i.rsplit('_',1)[0] + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                ts.legend.add_face(RectFace(width_dict[i], 20, '#FFFFFF', '#FFFFFF'), column=leg_column)
                for num2 in range(11):
                    y = num2 * 20 + 30
                    val = (max_val_dict[i] - min_val_dict[i]) * num2 / 10.0
                    h = val / (max_val_dict[i] - min_val_dict[i]) * 270
                    s = 0.5
                    l = 0.5
                    colour = hsl_to_str(h, s, l)
                    ts.legend.add_face(TextFace(font_gap * ' ' + str(val) + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                    ts.legend.add_face(RectFace(width_dict[i], 20, colour, colour), column=leg_column)
                leg_column += 2
            elif type_dict[i] == 'colour_scale_date':
                ts.legend.add_face(TextFace(font_gap * ' ' + i.rsplit('_',1)[0] + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                ts.legend.add_face(RectFace(width_dict[i], 20, '#FFFFFF', '#FFFFFF'), column=leg_column)
                for num2 in range(11):
                    y = num2 * 20 + 30
                    val = (max_val_dict[i] - min_val_dict[i]) * num2 / 10.0
                    h = val / (max_val_dict[i] - min_val_dict[i]) * 360
                    s = 0.5
                    l = 0.5
                    colour = hsl_to_str(h, s, l)
                    days = str(int(val / 60/ 60/ 24)) + ' days'
                    ts.legend.add_face(TextFace(font_gap * ' ' + days + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=leg_column+1)
                    ts.legend.add_face(RectFace(width_dict[i], 20, colour, colour), column=leg_column)
                leg_column += 2
            for n in t.traverse():
                if n.is_leaf():
                    name = leaf_name_dict[n.name]
                    if i in val_dict[name]:
                        val = val_dict[name][i]
                    else:
                        val = 'empty'
                    if type_dict[i] == 'colour':
                        n.add_face(RectFace(width_dict[i], 20, colour_dict[i][val], colour_dict[i][val]), column=num+1, position="aligned")
                    elif type_dict[i] == 'colour_scale' or type_dict[i] == 'colour_scale_date':
                        if val == 'empty':
                            colour = '#FFFFFF'
                        else:
                            h = (val - min_val_dict[i]) / (max_val_dict[i] - min_val_dict[i]) * 360
                            s = 0.5
                            l = 0.5
                            colour = hsl_to_str(h, s, l)
                        n.add_face(RectFace(width_dict[i], 20, colour, colour), column=num+1, position="aligned")
                    elif type_dict[i] == 'text':
                        n.add_face(TextFace(font_gap * ' ' + val + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=num+1, position="aligned")
    if not pres_abs is None:
        starting_col = len(column_list) + 1
        subprocess.Popen('makeblastdb -out tempdb -dbtype prot -in ' + pres_abs[0], shell=True).wait()
        folder = pres_abs[1]
        len_dict = {}
        gene_list = []
        ts.legend.add_face(TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=starting_col+1)
        ts.legend.add_face(RectFace(20, 20, '#FFFFFF', '#FFFFFF'), column=starting_col)
        ts.legend.add_face(TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=starting_col+1)
        ts.legend.add_face(RectFace(20, 20, "#5ba965", "#5ba965"), column=starting_col)
        ts.legend.add_face(TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True), column=starting_col+1)
        ts.legend.add_face(RectFace(20, 20, "#cb5b4c", "#cb5b4c"), column=starting_col)
        with open(pres_abs[0]) as f:
            for line in f:
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    gene_list.append(name)
                    len_dict[name] = 0
                    nameF = TextFace(font_gap * ' ' + name + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True)
                    nameF.rotation = -90
                    ts.aligned_header.add_face(nameF, column = starting_col +len(gene_list)-1)
                else:
                    len_dict[name] += len(line.rstrip())
        min_length = 0.9
        min_ident = 90
        for n in t.iter_leaves():
            the_name = n.name
            if the_name[0] == '"' and the_name[-1] == '"':
                the_name = the_name[1:-1]
            if the_name.endswith('.ref'):
                the_name = the_name[:-4]
            if not os.path.exists(folder + '/' + the_name):
                for q in os.listdir(folder):
                    if q.startswith(the_name):
                        the_name = q
            if not os.path.exists(the_name + '.blast'):
                subprocess.Popen('blastx -query ' + folder + '/' + the_name + ' -db tempdb -outfmt 6 -num_threads 24 -out ' + the_name + '.blast', shell=True).wait()
            gotit = set()
            with open(the_name + '.blast') as b:
                for line in b:
                    query, subject, ident, length = line.split()[:4]
                    ident = float(ident)
                    length = int(length)
                    if ident >= min_ident and length >= min_length * len_dict[subject]:
                        gotit.add(subject)
            for num, i in enumerate(gene_list):
                if i in gotit:
                    colour = "#5ba965"
                else:
                    colour = "#cb5b4c"
                n.add_face(RectFace(20, 20, colour, colour), column=num+starting_col, position="aligned")
        # for num, i in enumerate(gene_list):
        #     x = (starting_col + num) * 200
        #     svg.writeString(i, x+50, 20, 12)
        #     y = 30
        #     svg.drawOutRect(x + 50, y, 12, 12, strtorgb('#5ba965'), strtorgb('#5ba965'), lt=0)
        #     svg.writeString('present', x + 70, y + 12, 12)
        #     y = 50
        #     svg.drawOutRect(x + 50, y, 12, 12, strtorgb('#cb5b4c'), strtorgb('#cb5b4c'), lt=0)
        #     svg.writeString('absent', x + 70, y + 12, 12)




    # Set these to False if you don't want bootstrap/distance values
    ts.show_branch_length = label
    ts.show_branch_support = bootstrap
    ts.show_leaf_name = False
    for node in t.traverse():
        if node.is_leaf():
            node.add_face(AttrFace("name", fsize=font_size, ftype=font_type, tight_text=True, fgcolor='black'), column=0, position="aligned")

    ts.margin_left = 20
    ts.margin_right = 100
    ts.margin_top = 20
    ts.margin_bottom = 20
    if extend:
        ts.draw_guiding_lines = True
    ts.scale = the_scale
    if not circular is None:
        ts.mode = "c"
        ts.arc_start = 0
        ts.arc_span = 360
    if out_file is None:
        t.show(tree_style=ts)
    else:
        t.render(out_file, w=210, units='mm', tree_style=ts)


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--newick", help="newick tree file", required=True, metavar="tree.nw")
parser.add_argument("-o", "--output_svg", help="SVG of tree", metavar="tree.svg")
parser.add_argument("-c", '--color', action="store_true", default=False, help="Color tree")
parser.add_argument("-z", '--blend', action="store_true", default=False, help="Blend parent nodes based on child node colour.")
parser.add_argument("-l", '--label', action="store_true", default=False, help="Add distance values")
parser.add_argument("-b", '--bootstrap', action="store_true", default=False, help="add bootstrap values")
parser.add_argument("-s", "--scale", type=float, default=5000, help="x scale of tree")
parser.add_argument("-d", '--desc', help="Which columns and how to draw tree labels", metavar="assemblies.csv")
parser.add_argument("-t", '--tsv', help="add resistances from pathogendb", metavar="resistance.csv")
parser.add_argument("-g", '--group_file', help="file with groups of strains", metavar="Group.tsv")
parser.add_argument("-e", '--extend', help="extend tree branch", default=False, action="store_true")
parser.add_argument("-p", '--presence_absence', help="tblastn a fasta of proteins against each strain and report presence\
 or absense, give a protein fasta file, the folder in which the fastas used to generate the tree are located and a working directory.",
                    nargs=3, metavar=("genes.faa", "fasta_folder", "working_dir"))
parser.add_argument("-a", '--arc', type=int, help="To make a circular tree input arc angle. (i.e. 180, 360)")
args = parser.parse_args()

draw_tree(args.newick, args.color, args.blend, args.label, args.output_svg, args.scale, args.extend, args.bootstrap, args.group_file, args.desc, args.tsv, args.presence_absence, args.arc)
