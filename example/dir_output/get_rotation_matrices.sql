
create temp view tv_group_0 as
    select ball_id_DB, ball_id_Query, max(score) as max_score
    from match_rankings
    group by ball_id_DB, ball_id_Query;
 
create temp view tv_group_1 as
    select  
    match_rankings.PDB_pair as PDB_pair,
    match_rankings.ball_id_DB as ball_id_DB,
        match_rankings.ball_id_Query as ball_id_Queryball_id_Query,
        match_rankings.score as score,
		match_rankings.T0_0  as T0_0,
		match_rankings.T0_1  as T0_1,
		match_rankings.T0_2  as T0_2,
		match_rankings.R_00  as R_00,
		match_rankings.R_01  as R_01,
		match_rankings.R_02  as R_02,
		match_rankings.R_10  as R_10,
		match_rankings.R_11  as R_11,
		match_rankings.R_12  as R_12,
		match_rankings.R_20  as R_20,
		match_rankings.R_21  as R_21,
		match_rankings.R_22  as R_22,
		match_rankings.T1_0  as T1_0,
		match_rankings.T1_1  as T1_1,
		match_rankings.T1_2  as T1_2
    from  match_rankings, tv_group_0
    where match_rankings.ball_id_DB = tv_group_0.ball_id_DB AND
        match_rankings.ball_id_Query = tv_group_0.ball_id_Query AND
        match_rankings.score = tv_group_0.max_score;

select  T0_0,
        T0_1,
        T0_2,
        R_00,
        R_01,
        R_02,
        R_10,
        R_11,
        R_12,
        R_20,
        R_21,
        R_22,
        T1_0,
        T1_1,
        T1_2
from tv_group_1           
order by (0 - score);
