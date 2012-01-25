function term_reason = Two_method_crossover(control_mask, first_method_term_cond, first_search_method, second_method_term_cond, second_search_method)

first_search_method(control_mask, first_method_term_cond);

term_reason = second_search_method(control_mask, second_method_term_cond);