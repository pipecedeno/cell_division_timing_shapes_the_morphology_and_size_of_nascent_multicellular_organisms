
#Date: 5 may 2025
# This code is to add some missing columns that are needed to execute the synchrony_mother_daughter_31jan2024.csv
# code

df_all=read.csv("~/Desktop/timelapses_article/timelapse_joined_files_2025may5.csv", header = TRUE)

df_all$minutes=df_all$DELTA_T/60
df_all$hours=df_all$minutes/60
df_all$mother_delta_hours=df_all$mother_delta_t/3600

df_all$division_number=as.factor(df_all$division_number+1)

df_all$transfer=as.numeric(str_replace(df_all$timepoint,'t',''))

df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

#adding variable to make a unique id for each branch
df_all$unique_id_branch=paste(df_all$complete_id,df_all$TRACK_ID,sep='_')

write.csv(df_all, file="~/Desktop/timelapses_article/timelapse_doubling_time_inf_2025may5.csv", row.names = FALSE)
