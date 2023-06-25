use serde::{Serialize, Deserialize};
use sqlx::types::time::{PrimitiveDateTime};
use uuid::Uuid;

pub mod pb {
    tonic::include_proto!("auth");
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Account {
    pub id: String,
    pub user_id: Uuid,
    pub provider_type: String,
    pub provider_id: String,
    pub provider_account_id: String,
    pub refresh_token: Option<String>,
    pub access_token: Option<String>,
    pub access_token_expires: Option<PrimitiveDateTime>,
    pub created_at: PrimitiveDateTime,
    pub updated_at: PrimitiveDateTime,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Session {
    pub id: String,
    pub user_id: Uuid,
    pub expires: PrimitiveDateTime,
    pub session_token: String,
    pub access_token: String,
    pub created_at: PrimitiveDateTime,
    pub updated_at: PrimitiveDateTime,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct User {
    pub id: Uuid,
    pub name: Option<String>,
    pub email: Option<String>,
    pub email_verified: Option<PrimitiveDateTime>,
    pub image: Option<String>,
    pub created_at: PrimitiveDateTime,
    pub updated_at: PrimitiveDateTime,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VerificationRequest {
    pub id: String,
    pub identifier: String,
    pub token: String,
    pub expires: PrimitiveDateTime,
    pub created_at: PrimitiveDateTime,
    pub updated_at: PrimitiveDateTime,
}

fn main() {
    println!("Hello, world!");
}
