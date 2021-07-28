from azure.storage.blob import BlockBlobService, PublicAccess
blob_service_client = BlockBlobService(account_name='datasetsnpeff', sas_token='sv=2019-10-10&st=2020-09-01T00%3A00%3A00Z&se=2050-09-01T00%3A00%3A00Z&si=prod&sr=c&sig=isafOa9tGnYBAvsXFUMDGMTbsG2z%2FShaihzp7JE5dHw%3D')

# blob_service_client.list_blob_names("dataset") # a generator object

blob_service_client.get_blob_to_path('dataset/v5_0', "snpEff_v5_0_GRCh38.99.zip", "snpEff_v5_0_GRCh38.99.zip")

print("snpEff_v5_0_GRCh38.99.zip has been successfully downloaded.")
