classdef aFMM_Tree < handle
    % aFMM_Tree Adaptive FMM tree. We assume that the source and target points are all in the square
    % [0, 1]^2.
    
    % Reference:
    % J. Carrier, L. Greengard, and V. Rokhlin. A fast adaptive multipole algorithm for particle 
    % simulations. SIAM journal on scientific and statistical computing, 9(4):669–686, 1988.
    
    % Jingyu Liu, November 17, 2022.
    
    % TODO: Find why aFMM is slow than uniformFMM and fix it.
    
    properties
        % Tree information.
        level_ = 0;
        order_ = 0;
        num_level_ = 0;  % In fact, it is the highest level in the subtree whose root is the box.
        leaf_ = 0;  % Whether it is a leaf box.
        parent_;
        children_;
        neighbor_list_ = {};  % Adjacent boxes in the same level.
        U_list_ = {};  % The box itself and leaf boxes that are adjacent to the box when the box is a leaf, otherwise empty.
        V_list_ = {};  % Children of the neighbors of the box's parent which are well-separated from the box.
        W_list_ = {};  % Descendants of the box's neighbors whose parents are adjacent to the box but who are not adjacent to the box when the box is a leaf, otherwise empty. 
        X_list_ = {};  % Boxes whose W-list containes the box. 
        % Box and point information.
        center_ = zeros(1, 2);
        half_length_ = 0;
        source_points_ = [];
        source_charges_ = [];
        source_order_ = [];
        target_points_ = [];
        target_order_ = [];
        % FMM information.
        multipole_expansion_ = [];
        local_expansion_ = [];
        parent_local_expansion_ = [];
        
    end
    
    methods
        % Initialization.
        function obj = aFMM_Tree(source_points, source_charges, ...
                min_points, ...
                source_order, ...
                level, order, ...
                center, half_length)
            % aFMM_TREE Constructor.
            arguments
                source_points(:, 2);
                source_charges(:, 1);
                min_points(1, 1) = 64;
                source_order(:, 1) = 1 : length(source_charges);
                level(1, 1) = 0;
                order(1, 1) = 0;
                % The box must be square.
                center(1, 2) = [0.5, 0.5];
                half_length(1, 1) = 0.5;
            end
            
            obj.source_points_ = source_points;
            obj.source_charges_ = source_charges;
            obj.source_order_ = source_order;
            obj.level_ = level;
            obj.order_ = order;
            obj.center_ = center;
            obj.half_length_ = half_length;
            if obj.level_ == 0
                BuildTree(obj, min_points);
                SetList(obj);
            end
            
        end
        
        function BuildTree(obj, min_points)
            % BuilTree Build an FMM tree.
            
            n = length(obj.source_charges_);
            if n <= min_points
                obj.num_level_ = obj.level_;
                obj.leaf_ = 1;
                return;
            end
            
            % Partition the box into 4 sub-boxes.
            half_length = obj.half_length_ / 2;
            dxdy = [-half_length, half_length;
                half_length, half_length;
                -half_length, -half_length;
                half_length, -half_length];
            total_source_number = 0;
            for iter = 1 : 4
                center = obj.center_ + dxdy(iter, :);
                
                source_index = 1 : n;
                source_index = intersect(source_index, ...
                    find(obj.source_points_(:, 1) >= center(1, 1) - half_length));
                source_index = intersect(source_index, ...
                    find(obj.source_points_(:, 1) < center(1, 1) + half_length));
                source_index = intersect(source_index, ...
                    find(obj.source_points_(:, 2) >= center(1, 2) - half_length));
                source_index = intersect(source_index, ...
                    find(obj.source_points_(:, 2) < center(1, 2) + half_length));
                source_points = obj.source_points_(source_index, :);
                source_charges = obj.source_charges_(source_index);
                source_order = obj.source_order_(source_index);
                total_source_number = total_source_number + length(source_index);
                
                obj.children_{iter} = aFMM_Tree(source_points, source_charges, ...
                    min_points, ...
                    source_order, ...
                    obj.level_ + 1, 4 * obj.order_ + iter - 1, ...
                    center, half_length);
                obj.children_{iter}.parent_ = obj;
            end
            
            if total_source_number ~= n
                error("Source points appear on the boundary!")
            end
            
            % Recursively buildtree.
            for iter = 1 : 4
                BuildTree(obj.children_{iter}, min_points);
            end
            
            obj.num_level_ = max(...
                [obj.children_{1}.num_level_, ...
                obj.children_{2}.num_level_, ...
                obj.children_{3}.num_level_, ...
                obj.children_{4}.num_level_]);
            
        end
        
        function SetList(obj)
            % SetList Set the FMM lists.
            
            if obj.leaf_ == 1
                % Update U-list.
                parent = obj.parent_;
                if isempty(parent)
                    return;
                end
                for iter = 1 : 4
                    UpdateUList(obj, parent.children_{iter});
                end
                for i = 1 : length(parent.neighbor_list_)
                    UpdateUList(obj, parent.neighbor_list_{i});
                end
                
                % Update W-list and X-list.
                for i = 1 : length(obj.neighbor_list_)
                    nb_box = obj.neighbor_list_{i};
                    if nb_box.leaf_ == 0
                        for t = 1 : 4
                            UpdateWXList(obj, nb_box.children_{t});
                        end
                    end
                end
                
                return;
            end
            
            % Update neighbor list and V-list.
            for iter = 1 : 4
                child = obj.children_{iter};
                for it = 1 : iter - 1
                    child.neighbor_list_{end + 1} = obj.children_{it};
                end
                for it = iter + 1 : 4
                    child.neighbor_list_{end + 1} = obj.children_{it};
                end
                for i = 1 : length(obj.neighbor_list_)
                    nb_box = obj.neighbor_list_{i};
                    if nb_box.leaf_ == 0
                        for nb_iter = 1 : 4
                            nb_box_child = nb_box.children_{nb_iter};
                            % nb_box_child and child are of the same level.
                            if max(abs(nb_box_child.center_ - child.center_)) <= obj.half_length_
                                % obj.half_length_ = child.half_length_ + nb_box_child.half_length_.
                                child.neighbor_list_{end + 1} = nb_box_child;
                            else
                                child.V_list_{end + 1} = nb_box_child;
                            end
                        end
                    end
                end
            end
            
            % Recursively.
            for iter = 1 : 4
                SetList(obj.children_{iter});
            end
            
        end
        
        function UpdateUList(obj, box)
            % UpdateUList Update U-list of obj.
                       
            if box.leaf_ == 1
                if obj.level_ > box.level_ || ...
                        (obj.level_ == box.level_ && obj.order_ >= box.order_)
                    % Naively, one may directly set the boxes in U-list, but there is a problem: The
                    % smaller box can't find the larger box in its U-list. Therefore, we use the
                    % following property: If B is in the U-list of A, then A is also in the U-list 
                    % of B.
                    return;
                end
                if max(abs(box.center_ - obj.center_)) <= obj.half_length_ + box.half_length_
                    obj.U_list_{end + 1} = box;
                    box.U_list_{end + 1} = obj;
                end
            elseif max(abs(box.center_ - obj.center_)) <= obj.half_length_ + box.half_length_
                % If a leaf box is adjacent to obj, its parent must be adjacent to obj.
                for iter = 1 : 4
                    UpdateUList(obj, box.children_{iter});
                end
            end
            
        end
        
        function UpdateWXList(obj, box)
            % UpdateWXList Update W-list and X-list of obj.
            
            box_parent = box.parent_;  % It must be nonempty!
            if max(abs(box.center_ - obj.center_)) > obj.half_length_ + box.half_length_
                % box's children can't be in the W-list of obj.
                if max(abs(box_parent.center_ - obj.center_)) <= obj.half_length_ + box_parent.half_length_
                    % box_parent is adjacent to obj and box is not adjacent to obj.
                    obj.W_list_{end + 1} = box;
                    box.X_list_{end + 1} = obj;
                    return;
                end
            else
                if box.leaf_ == 0
                    for iter = 1 : 4
                        UpdateWXList(obj, box.children_{iter});
                    end
                end
            end
            
        end
        
        % FMM algorithm.
        function FMM_Alg(obj, tol)
            % FMM_Alg FMM algorithms.
            p = ceil(-log2(tol));  % Number of expansion terms.
            
            global M2M_combination;
            M2M_combination = zeros(p, p);  % (l - 1, k - 1) for 1 <= k <= l <= p. 
            for l = 1 : p
                for k = 1 : l
                    M2M_combination(l, k) = nchoosek(l - 1, k - 1);
                end
            end
            global M2L_combination;
            M2L_combination = zeros(p, p);  % (l + k - 1, k - 1) for 1 <= k, l <= p.
            for l = 1 : p
                for k = 1 : p
                    M2L_combination(l, k) = nchoosek(l + k - 1, k - 1);
                end
            end
            global L2L_combination;
            L2L_combination = zeros(p + 1, p + 1);  % (k, l) for 0 <= l <= k <= p.
            for l = 0 : p
                for k = l : p
                    L2L_combination(k + 1, l + 1) = nchoosek(k, l);
                end
            end
            
            % S2M and M2M.
            for tmplevel = obj.num_level_ : -1 : 0
                RecursiveS2M2M(obj, p, tmplevel);
            end
            
            % Set parent local expansion zero for boxes in level 1.
            if obj.leaf_ == 0
                for iter = 1 : 4
                    obj.children_{iter}.parent_local_expansion_ = zeros(p + 1, 1);
                end
            end
            
            % M2L and L2L.
            for tmplevel = 1 : obj.num_level_
                RecursiveM2L2L(obj, p, tmplevel);
            end
            
        end
        
        function RecursiveS2M2M(obj, p, whatlevel)
            % RecursiveS2M2M Source to multipole and multipole to multipole recursively.
            
            if obj.level_ == whatlevel
                S2M2M(obj, p);
            else
                if obj.leaf_ == 0
                    for iter = 1 : 4
                        RecursiveS2M2M(obj.children_{iter}, p, whatlevel);
                    end
                end
            end
            
        end
        
        function S2M2M(obj, p)
            % S2M2M Source to multipole and multipole to multipole.
            
            if obj.leaf_ == 1
                S2M(obj, p);
            else
                M2M(obj, p);
            end
            
        end
        
        function S2M(obj, p)
            % S2M Source to multipole.
            
            m = length(obj.source_charges_);
            if m == 0
                obj.multipole_expansion_ = zeros(p + 1, 1);
                return;
            end
            z = complex(obj.source_points_(:, 1), obj.source_points_(:, 2));
            zc = complex(obj.center_(1, 1), obj.center_(1, 2));
            obj.multipole_expansion_ = zeros(p + 1, 1);
            obj.multipole_expansion_(1) = sum(obj.source_charges_);
%             for k = 1 : p
%                 for i = 1 : m
%                     obj.multipole_expansion_(k + 1) = obj.multipole_expansion_(k + 1) + ...
%                         (-obj.source_charges_(i) * (z(i) - zc)^k);
%                 end
%                 obj.multipole_expansion_(k + 1) = obj.multipole_expansion_(k + 1) / k;
%             end
            z_minus_zc = z - zc;
            zexp = z_minus_zc;
            for k = 1 : p
                obj.multipole_expansion_(k + 1) = -zexp.' * obj.source_charges_ / k;
                zexp = zexp .* z_minus_zc;
            end
            
        end
        
        function M2M(obj, p)
            % M2M Multipole to multipole.
            
            global M2M_combination;
            
            obj.multipole_expansion_ = zeros(p + 1, 1);
            zc = complex(obj.center_(1, 1), obj.center_(1, 2));
            for iter = 1 : 4
                child = obj.children_{iter};
                z0 = complex(child.center_(1, 1), child.center_(1, 2));
                obj.multipole_expansion_(1) = obj.multipole_expansion_(1) + child.multipole_expansion_(1);
                for l = 1 : p
                    for k = 1 : l
                        obj.multipole_expansion_(l + 1) = obj.multipole_expansion_(l + 1) + ...
                            child.multipole_expansion_(k + 1) * (z0 - zc)^(l - k) * M2M_combination(l, k);
                    end
                    obj.multipole_expansion_(l + 1) = obj.multipole_expansion_(l + 1) - ...
                        child.multipole_expansion_(1) * (z0 - zc)^l / l;
                end
            end
            
        end
        
        function RecursiveM2L2L(obj, p, whatlevel)
            % RecursiveM2L2L Multipole to local and local to local recursively.
            
            if obj.level_ == whatlevel
                M2L2L(obj, p);
            else
                if obj.leaf_ == 0
                    for iter = 1 : 4
                        RecursiveM2L2L(obj.children_{iter}, p, whatlevel);
                    end
                end
            end
        end
        
        function M2L2L(obj, p)
            % M2L2L Multipole to local and local to local.
            
            M2L(obj, p);
            
            if obj.leaf_ == 0
                L2L(obj, p);
            end
            
        end
        
        function M2L(obj, p)
            % M2L Multipole to local.
            
            global M2L_combination;
            
            obj.local_expansion_ = zeros(p + 1, 1);
            zc = complex(obj.center_(1, 1), obj.center_(1, 2));
            % M2L from V-list.
            for i = 1 : length(obj.V_list_)
                V_box = obj.V_list_{i};
                z0 = complex(V_box.center_(1, 1), V_box.center_(1, 2));
                obj.local_expansion_(1) = obj.local_expansion_(1) + ...
                    V_box.multipole_expansion_(1) * log(-(z0 - zc));
                for l = 1 : p
                    obj.local_expansion_(l + 1) = obj.local_expansion_(l + 1) + ...
                        (-V_box.multipole_expansion_(1) / l / (z0 - zc)^l);
                end
                for k = 1 : p
                    obj.local_expansion_(1) = obj.local_expansion_(1) + ...
                        V_box.multipole_expansion_(k + 1) / (z0 - zc)^k * (-1)^k;
                    for l = 1 : p
                        obj.local_expansion_(l + 1) = obj.local_expansion_(l + 1) + ...
                            V_box.multipole_expansion_(k + 1) / (z0 - zc)^k * M2L_combination(l, k) * (-1)^k / (z0 - zc)^l;
                    end
                end
            end
            % M2L from X-list.
            for i = 1 : length(obj.X_list_)
                X_box = obj.X_list_{i};
                m = length(X_box.source_charges_);
                z0 = complex(X_box.source_points_(:, 1), X_box.source_points_(:, 2));
                for t = 1 : m
                    obj.local_expansion_(1) = obj.local_expansion_(1) + ...
                        X_box.source_charges_(t) * log(-(z0(t) - zc));
                    for l = 1 : p
                        obj.local_expansion_(l + 1) = obj.local_expansion_(l + 1) + ...
                            (-X_box.source_charges_(t) / l / (z0(t) - zc)^l);
                    end
                end
            end
            % Combine together.
            obj.local_expansion_ = obj.local_expansion_ + obj.parent_local_expansion_;
            
        end
        
        function L2L(obj, p)
            % L2L Local to local.
            
            global L2L_combination;
            
            z0 = complex(obj.center_(1, 1), obj.center_(1, 2));
            for iter = 1 : 4
                child = obj.children_{iter};
                child.parent_local_expansion_ = zeros(p + 1, 1);
                zc = complex(child.center_(1, 1), child.center_(1, 2));
                for l = 0 : p
                    for k = l : p
                        child.parent_local_expansion_(l + 1) = child.parent_local_expansion_(l + 1) + ...
                            obj.local_expansion_(k + 1) * L2L_combination(k + 1, l + 1) * (-(z0 - zc))^(k - l);
                    end
                end
            end
            
        end
        
        % FMM computation.
        function potential = FMM_Compute(obj, target_points)
            % FMM_Compute Compute the potentials on the target points using FMM.
            
            m = size(target_points, 1);
            potential = zeros(m, 1);
            obj.target_points_ = target_points;
            obj.target_order_ = 1 : m;
            potential = RecursiveL2T(obj, potential);
            potential = real(-potential);
            
        end
        
        function potential = RecursiveL2T(obj, potential)
            % RecursiveL2T Local to target recursively.
            
            if obj.leaf_ == 1
                potential = L2T(obj, potential);
                % Clear.
                obj.target_points_= [];
                obj.target_order_ = [];
                return;
            end
            
            m = size(obj.target_points_, 1);
            total_target_number = 0;
            for iter = 1 : 4
                child = obj.children_{iter};
                center = child.center_;
                half_length = child.half_length_;
                
                target_index = 1 : m;
                target_index = intersect(target_index, ...
                    find(obj.target_points_(:, 1) >= center(1, 1) - half_length));
                target_index = intersect(target_index, ...
                    find(obj.target_points_(:, 1) < center(1, 1) + half_length));
                target_index = intersect(target_index, ...
                    find(obj.target_points_(:, 2) >= center(1, 2) - half_length));
                target_index = intersect(target_index, ...
                    find(obj.target_points_(:, 2) < center(1, 2) + half_length));
                target_points = obj.target_points_(target_index, :);
                target_order = obj.target_order_(target_index);
                total_target_number = total_target_number + length(target_index);
                
                child.target_points_ = target_points;
                child.target_order_ = target_order;
            end
            
            if total_target_number ~= m
                error("Target points appear on the boundary!")
            end
            
            % Recursively.
            for iter = 1 : 4
                potential = RecursiveL2T(obj.children_{iter}, potential);
            end
            
            % Clear.
            obj.target_points_= [];
            obj.target_order_ = [];
            
        end
        
        function potential = L2T(obj, potential)
            % Local to target.
            
            z = complex(obj.target_points_(:, 1), obj.target_points_(:, 2));
            % Far-field.
            z0 = complex(obj.center_(1, 1), obj.center_(1, 2));
            potential(obj.target_order_) = potential(obj.target_order_) + ...
                polyval(flip(obj.local_expansion_), z - z0);
            % Near-field: The box itself and its U-list.
            zz = complex(obj.source_points_(:, 1), obj.source_points_(:, 2));
%             for k = 1 : length(obj.source_charges_)
%                 zz = complex(obj.source_points_(k, 1), obj.source_points_(k, 2));
%                 potential(obj.target_order_) = potential(obj.target_order_) + ...
%                     obj.source_charges_(k) * log(z - zz);
%             end
            potential(obj.target_order_) = potential(obj.target_order_) + ...
                sum(obj.source_charges_' .* log(kron(z, ones(size(zz.'))) - kron(ones(size(z)), zz.')), 2);
            for i = 1 : length(obj.U_list_)
                U_box = obj.U_list_{i};
                zz = complex(U_box.source_points_(:, 1), U_box.source_points_(:, 2));
%                 for k = 1 : length(U_box.source_charges_)
%                     zz = complex(U_box.source_points_(k, 1), U_box.source_points_(k, 2));
%                     potential(obj.target_order_) = potential(obj.target_order_) + ...
%                         U_box.source_charges_(k) * log(z - zz);
%                 end
                potential(obj.target_order_) = potential(obj.target_order_) + ...
                    sum(U_box.source_charges_' .* log(kron(z, ones(size(zz.'))) - kron(ones(size(z)), zz.')), 2);
            end
            % Medium-field: W-list.
            for i = 1 : length(obj.W_list_)
                W_box = obj.W_list_{i};
                z0 = complex(W_box.center_(1, 1), W_box.center_(1, 2));
                potential(obj.target_order_) = potential(obj.target_order_) + ...
                    W_box.multipole_expansion_(1) * log(z - z0) + ...
                    polyval([flip(W_box.multipole_expansion_(2 : end)); 0], 1 ./ (z - z0));
            end
            
        end
        
    end
    
end