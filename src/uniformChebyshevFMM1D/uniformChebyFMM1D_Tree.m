classdef uniformChebyFMM1D_Tree < handle
    % uniformChebyFMM1D_Tree. 1D Uniform Chebyshev FMM tree. We assume that the source and target 
    % points are all in the interval [0, 1].
    
    % Reference:
    % W. Fong and E. Darve. The black-box fast multipole method. Journal of computational physics, 
    % 228(23):8712–8725, 2009.
    
    % Jingyu Liu, November 24, 2022.
    
    properties
        % Tree information.
        level_ = 0;
        order_ = 0;
        num_level_ = 0;  % In fact, it is the highest level in the subtree whose root is the box.
        leaf_ = 0;  % Whether it is a leaf box.
        parent_;
        children_;
        neighbor_list_ = {};  % Adjacent boxes in the same level.
        interaction_list_ = {};  % Children of the neighbors of the box's parent which are well-separated from the box.
        % Box and point information.
        center_ = 0.5;
        half_length_ = 0.5;
        source_points_ = [];
        source_charges_ = [];
        source_order_ = [];
        cheby_points_ = [];  % Chebyshev points in the interval.
        target_points_ = [];
        target_order_ = [];
        % FMM information.
        multipole_expansion_ = [];
        local_expansion_ = [];
        parent_local_expansion_ = [];
        
    end
    
    methods
        % Initialization.
        function obj = uniformChebyFMM1D_Tree(source_points, source_charges, source_order, ...
                level, order, ...
                center, half_length)
            % uniformChebyFMM1D_Tree Constructor.
            arguments
                source_points(:, 1);
                source_charges(:, 1);
                source_order(:, 1) = 1 : length(source_charges);
                level(1, 1) = 0;
                order(1, 1) = 0;
                % The box must be square.
                center(1, 1) = 0.5;
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
                min_points = 512;  % The number of particles in each leaf box is nearly min_points when uniform case.
                num_level = max(ceil(log2(length(obj.source_charges_) / min_points) / 2), 0);
                BuildTree(obj, num_level);
                SetList(obj);
            end
            
        end
        
        function BuildTree(obj, num_level)
            % BuilTree Build an FMM tree.
            
            if obj.level_ == num_level
                obj.num_level_ = obj.level_;
                obj.leaf_ = 1;
                return;
            end
            
            % Partition the box into 2 sub-boxes.
            half_length = obj.half_length_ / 2;
            dx = [-half_length; half_length];
            total_source_number = 0;
            
            for iter = 1 : 2
                center = obj.center_ + dx(iter);
                
                source_index = 1 : length(obj.source_charges_);
                source_index = intersect(source_index, ...
                    find(obj.source_points_ >= center - half_length));
                source_index = intersect(source_index, ...
                    find(obj.source_points_ < center + half_length));
                source_points = obj.source_points_(source_index);
                source_charges = obj.source_charges_(source_index);
                source_order = obj.source_order_(source_index);
                total_source_number = total_source_number + length(source_index);
                
                obj.children_{iter} = uniformChebyFMM1D_Tree(source_points, source_charges, source_order, ...
                    obj.level_ + 1, 2 * obj.order_ + iter - 1, ...
                    center, half_length);
                obj.children_{iter}.parent_ = obj;
            end
            
            if total_source_number ~= length(obj.source_charges_)
                error("Source points appear on the boundary!")
            end
            
            % Recursively buildtree.
            for iter = 1 : 2
                BuildTree(obj.children_{iter}, num_level);
            end
            
            obj.num_level_ = max(obj.children_{1}.num_level_, obj.children_{2}.num_level_);
            
        end
        
        function SetList(obj)
            % SetList Set the FMM lists.
            
            if obj.leaf_ == 1
                % We stand on the parent box.
                return;
            end
            
            % Update neighbor list and interaction list.
            for iter = 1 : 2
                child = obj.children_{iter};
                child.neighbor_list_{end + 1} = obj.children_{3 - iter};
                for i = 1 : length(obj.neighbor_list_)
                    nb_box = obj.neighbor_list_{i};
                    if nb_box.leaf_ == 0
                        for nb_iter = 1 : 2
                            nb_box_child = nb_box.children_{nb_iter};
                            % nb_box_child and child are of the same level.
                            if abs(nb_box_child.center_ - child.center_) <= obj.half_length_
                                % obj.half_length_ = child.half_length_ + nb_box_child.half_length_.
                                child.neighbor_list_{end + 1} = nb_box_child;
                            else
                                child.interaction_list_{end + 1} = nb_box_child;
                            end
                        end
                    end
                end
            end
            
            % Recursively.
            for iter = 1 : 2
                SetList(obj.children_{iter});
            end
            
        end
        
        % FMM algorithm.
        function FMM_Alg(obj, kernel_fun, tol)
            % FMM_Alg FMM algorithms.
            
            p = ceil(-log2(tol));  % Number of Chebyshev points in each interval.
            
            global kernel_function;
            kernel_function = kernel_fun;
            
            % S2M and M2M.
            for tmplevel = obj.num_level_ : -1 : 0
                RecursiveS2M2M(obj, p, tmplevel);
            end
            
            % Set parent local expansion zero for boxes in level 1.
            if obj.leaf_ == 0
                for iter = 1 : 2
                    obj.children_{iter}.parent_local_expansion_ = zeros(p, 1);
                end
            else
                obj.local_expansion_ = zeros(p, 1);
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
                for iter = 1 : 2
                    RecursiveS2M2M(obj.children_{iter}, p, whatlevel);
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
            
            left_point = obj.center_ - obj.half_length_;
            right_point = obj.center_ + obj.half_length_;
            obj.cheby_points_ = ChebyPoints(left_point, right_point, p);
            
            m = length(obj.source_charges_);
            if m == 0
                obj.multipole_expansion_ = zeros(p, 1);
                return;
            end
            
%             obj.multipole_expansion_ = zeros(1, p);
%             for k = 1 : p
%                 for j = 1 : m
%                     obj.multipole_expansion_(k) = obj.multipole_expansion_(k) + ...
%                         obj.source_charges_(j) * ...
%                         ComputeSp(left_point, right_point, p, obj.cheby_points_(k), obj.source_points_(j));
%                 end
%             end
            obj.multipole_expansion_ = ...
                ComputeSp(left_point, right_point, p, obj.cheby_points_, obj.source_points_) * ...
                obj.source_charges_;
            
        end
        
        function M2M(obj, p)
            % M2M Multipole to multipole.
            
            left_point = obj.center_ - obj.half_length_;
            right_point = obj.center_ + obj.half_length_;
            obj.cheby_points_ = ChebyPoints(left_point, right_point, p);
            
            obj.multipole_expansion_ = zeros(p, 1);
            for iter = 1 : 2
                child = obj.children_{iter};
%                 for k = 1 : p
%                     for l = 1 : p
%                         obj.multipole_expansion_(k) = obj.multipole_expansion_(k) + ...
%                             child.multipole_expansion_(l) * ...
%                             ComputeSp(left_point, right_point, p, obj.cheby_points_(k), child.cheby_points_(l));
%                     end
%                 end
                obj.multipole_expansion_ = obj.multipole_expansion_ + ...
                    ComputeSp(left_point, right_point, p, obj.cheby_points_, child.cheby_points_) * ...
                    child.multipole_expansion_;
            end
            
        end
        
        function RecursiveM2L2L(obj, p, whatlevel)
            % RecursiveM2L2L Multipole to local and local to local recursively.
            
            if obj.level_ == whatlevel
                M2L2L(obj, p);
            else
                for iter = 1 : 2
                    RecursiveM2L2L(obj.children_{iter}, p, whatlevel);
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
            
            global kernel_function;
            
            obj.local_expansion_ = zeros(p, 1);
            % M2L from interaction-list.
            for i = 1 : length(obj.interaction_list_)
                interact_box = obj.interaction_list_{i};
%                 for k = 1 : p
%                     for l = 1 : p
%                         obj.local_expansion_(k) = obj.local_expansion_(k) + ...
%                             interact_box.multipole_expansion_(l) * ...
%                             kernel_function(obj.cheby_points_(k), interact_box.cheby_points_(l));
%                     end
%                 end
            obj.local_expansion_ = obj.local_expansion_ + ...
                kernel_function(obj.cheby_points_, interact_box.cheby_points_') * ...
                interact_box.multipole_expansion_;
                    
            end
            % Combine together.
            obj.local_expansion_ = obj.local_expansion_ + obj.parent_local_expansion_;
            
        end
        
        function L2L(obj, p)
            % L2L Local to local.
            
            left_point = obj.center_ - obj.half_length_;
            right_point = obj.center_ + obj.half_length_;
            
            for iter = 1 : 2
                child = obj.children_{iter};
                child.parent_local_expansion_ = zeros(p, 1);
%                 for k = 1 : p
%                     for l = 1 : p
%                         child.parent_local_expansion_(k) = child.parent_local_expansion_(k) + ...
%                             obj.local_expansion_(l) * ...
%                             ComputeSp(left_point, right_point, p, child.cheby_points_(k), obj.cheby_points_(l));
%                     end
%                 end
                child.parent_local_expansion_ = child.parent_local_expansion_ + ...
                    ComputeSp(left_point, right_point, p, child.cheby_points_, obj.cheby_points_) * ...
                    obj.local_expansion_;
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
            
        end
        
        function potential = RecursiveL2T(obj, potential)
            % RecursiveL2T Local to target recursively.
            
            if obj.leaf_ == 1
                potential = L2T(obj, potential);
                % Clear.
                % obj.target_points_= [];
                % obj.target_order_ = [];
                return;
            end
            
            m = size(obj.target_points_, 1);
            total_target_number = 0;
            for iter = 1 : 2
                child = obj.children_{iter};
                center = child.center_;
                half_length = child.half_length_;
                
                target_index = 1 : m;
                target_index = intersect(target_index, ...
                    find(obj.target_points_ >= center - half_length));
                target_index = intersect(target_index, ...
                    find(obj.target_points_ < center + half_length));
                target_points = obj.target_points_(target_index);
                target_order = obj.target_order_(target_index);
                total_target_number = total_target_number + length(target_index);
                
                child.target_points_ = target_points;
                child.target_order_ = target_order;
            end
            
            if total_target_number ~= m
                error("Target points appear on the boundary!")
            end
            
            % Recursively.
            for iter = 1 : 2
                potential = RecursiveL2T(obj.children_{iter}, potential);
            end
            
            % Clear.
            % obj.target_points_= [];
            % obj.target_order_ = [];
            
        end
        
        function potential = L2T(obj, potential)
            % Local to target.
            
            global kernel_function;
            
            left_point = obj.center_ - obj.half_length_;
            right_point = obj.center_ + obj.half_length_;
            
            p = length(obj.local_expansion_);
            
            % Far-field.
%             for i = 1 : length(obj.target_order_)
%                 for k = 1 : p
%                     potential(obj.target_order_(i)) = potential(obj.target_order_(i)) + ...
%                         obj.local_expansion_(k) * ...
%                         ComputeSp(left_point, right_point, p, obj.target_points_(i), obj.cheby_points_(k));
%                 end
%             end
            potential(obj.target_order_) = potential(obj.target_order_) + ...
                ComputeSp(left_point, right_point, p, obj.target_points_, obj.cheby_points_) * ...
                obj.local_expansion_;
            % Near-field.
%             for i = 1 : length(obj.target_order_)
%                 for k = 1 : length(obj.source_charges_)
%                     potential(obj.target_order_(i)) = potential(obj.target_order_(i)) + ...
%                         obj.source_charges_(k) * ...
%                         kernel_function(obj.target_points_(i), obj.source_points_(k));
%                 end
%             end
            potential(obj.target_order_) = potential(obj.target_order_) + ...
                kernel_function(obj.target_points_, obj.source_points_') * obj.source_charges_;
            for t = 1 : length(obj.neighbor_list_)
                nb_box = obj.neighbor_list_{t};
%                 for i = 1 : length(obj.target_order_)
%                     for k = 1 : length(nb_box.source_charges_)
%                         potential(obj.target_order_(i)) = potential(obj.target_order_(i)) + ...
%                             nb_box.source_charges_(k) * ...
%                             kernel_function(obj.target_points_(i), nb_box.source_points_(k));
%                     end
%                 end
                potential(obj.target_order_) = potential(obj.target_order_) + ...
                                kernel_function(obj.target_points_, nb_box.source_points_') * nb_box.source_charges_;
            end
            
        end
        
    end
    
end

