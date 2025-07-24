function [W] = MakeW(nx_grid,ny_grid,alpha)

%Create the W matrix. This is the discrete Laplacian matrix 
%(https://en.wikipedia.org/wiki/Laplacian_matrix).

%Note: I am not currently dividing by cellsize^2, so changing the cell size
%without changing alpha will change the length-scale of smoothing.

%Find the number of non-zero elements in W.
%The number of entries for each point is itself plus its number of
%horizontal and vertical neighbours: So 5 for an interior point, 4 for an
%edge point, and 3 for a corner point.
ngrid = nx_grid*ny_grid;
ncorner = 4;
nedge = 2*(nx_grid-2) + 2*(ny_grid-2); %Non-corner edge. The -2 terms are to exclude the corners.
ninterior = ngrid - nedge - ncorner;
count = ninterior * 5 + nedge * 4 + ncorner * 3; %Total number of non-zero entries in W.

%Initialize the arrays to hold the non-zero entries and their indices.
W_i=zeros([count,1]);
W_j=zeros([count,1]);
W_v=zeros([count,1]);
W_alpha_v=zeros([count,1]);

%Build the W matrix.
fprintf(1,'Building W vectors \n');
count=0;
for j=1:ny_grid
    for i=1:nx_grid
        row_index=i+(j-1)*nx_grid;
        weight=4;
        %Look left
        i1=i-1;
        j1=j;
        if i1<1
            weight=weight-1;
        else
            count=count+1;
            col_index=i1+(j1-1)*nx_grid;
            W_i(count)=row_index;
            W_j(count)=col_index;
            W_v(count)=-1;
        end
        %Look right
        i2=i+1;
        j2=j;
        if i2>nx_grid
            weight=weight-1;
        else
            count=count+1;
            col_index=i2+(j2-1)*nx_grid;
            W_i(count)=row_index;
            W_j(count)=col_index;
            W_v(count)=-1;

        end
        %Look down
        i3=i;
        j3=j-1;
        if j3<1
            weight=weight-1;
        else
            count=count+1;
            col_index=i3+(j3-1)*nx_grid;
            W_i(count)=row_index;
            W_j(count)=col_index;
            W_v(count)=-1;
        end
        %Look up
        i4=i;
        j4=j+1;
        if j4>ny_grid
            weight=weight-1;
        else
            count=count+1;
            col_index=i4+(j4-1)*nx_grid;
            W_i(count)=row_index;
            W_j(count)=col_index;
            W_v(count)=-1;
        end
        count=count+1;
        col_index=i+(j-1)*nx_grid;
        W_i(count)=row_index;
        W_j(count)=col_index;
        W_v(count)=weight;
    end
end

fprintf(1,'Assembling W \n');

%Multiply W_v by alpha where alpha controls the damping.
W_alpha_v(1:count)=W_v(1:count)*alpha;

%Build W_alpha
W=sparse(W_i,W_j,W_alpha_v,(nx_grid*ny_grid),(nx_grid*ny_grid));

fprintf(1,'W built\n');

end