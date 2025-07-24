function A = MakeA(S,delta_chi,grid_size,nn_grid,nx_grid,x1,y1)

%The number of rows in A is the number of non-base nodes, since all of
%these have an elevation that can be fit for. The number of columns is the
%number of grid cells at which we are fitting for ksn.

%The number of elements of delta_chi should be equal to the total number of
%nodes in S (not just the non-base nodes) so it can be indexed with S.ix
%and S.ixc, but values for base-level nodes will not be used.

n_nodes = length(S.x);
n_non_base = length(S.ix); %Number of non-base-level nodes in S.
A_i = cell(n_nodes,1);
A_j = cell(n_nodes,1);
A_v = cell(n_nodes,1);
i1=floor((S.x-x1)/grid_size)+1; %Calculate indices of the Ksn grid cells in which each stream node lies.
j1=floor((S.y-y1)/grid_size)+1;
col_index = i1+(j1-1)*nx_grid;
for irow = n_non_base:-1:1 %Go through the stack from the base up.
    id = S.ix(irow); %Index in S of the node this row predicts z at.
    irec = S.ixc(irow); %Index in S of the receiver of node id.
    nnzrec = length(A_v{irec});
    if (col_index(id) == col_index(irec)) && (nnzrec>0)
        %Node and its receiver are in the same ksn grid cell.
        nnzi = nnzrec; %Number of non-zero values in the row.
        A_v{id} = A_v{irec};
        A_j{id} = A_j{irec};
    else
        %Node is in a different ksn grid cell from its receiver.
        nnzi = nnzrec+1;
        A_v{id} = zeros(nnzi,1);
        A_j{id} = zeros(nnzi,1);
        A_v{id}(1:end-1) = A_v{irec};
        A_j{id}(1:end-1) = A_j{irec};
        A_j{id}(end) = col_index(id);
    end
    A_v{id}(end) = A_v{id}(end) + delta_chi(id); %Add the new delta_chi value to the appropriate column.
    A_i{id}=irow*ones(nnzi,1);
end
A_i = vertcat(A_i{:});
A_j = vertcat(A_j{:});
A_v = vertcat(A_v{:});
%Note: Entries in A_i, A_j, and A_v corresponding to base-level nodes will
%be empty, but these will simply disappear when vertcat is applied.

fprintf(1,'A vectors defined, assembling sparse matrix \n');

A=sparse(A_i,A_j,A_v,n_non_base,nn_grid);

end